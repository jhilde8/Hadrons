/*
 * A2AMesonFieldCast.cpp, part of Hadrons
 *
 * Converts a tree of double-precision meson-field HDF5 files to single
 * precision. Every <src>/<rel>.h5 is written to <dst>/<rel>.h5. The source
 * tree is opened read-only and never modified.
 *
 * Per file, the HDF5 tree is copied object by object:
 *   - every complex (re/im compound) double dataset is rewritten as ComplexF
 *     (Hdf5Type<ComplexF>, the type a ComplexF production build writes), with
 *     the source's dataset creation property list, so the chunk dims and the
 *     filters (Fletcher32) are reproduced;
 *   - groups, attributes and all other datasets are copied unchanged, so the
 *     metadata of every family (mf, emf, cmf, ...) is carried over without a
 *     per-family metadata class.
 * The data is read in bands of whole chunks, (chunk t) x (chunk rows) x (all
 * columns), so files with one timeslice (timeSliceIO) and files with all
 * timeslices both work with bounded memory. A complex dataset that is already
 * single precision is an error, as is a non-finite value or one outside the
 * float range.
 *
 * Integrity: the source is read through its Fletcher32 filter, so a corrupt
 * chunk fails the read. Element-wise comparison of the output against the
 * source is done separately by tests/Test_mf_cast_compare.
 *
 * Output: <dst>/<rel>.h5.partial is written, closed, then hard-linked to
 * <dst>/<rel>.h5 and unlinked; link() fails if the target exists, so a
 * finished file is never overwritten. A file whose target exists is skipped,
 * so a rerun after a failure continues where it stopped. A stale .partial
 * left by a killed run is replaced when its file is converted again. This
 * assumes nothing but this tool writes into dst.
 *
 * Parallelism: no MPI. The sorted file list is split round-robin by the
 * launcher's task rank (SLURM_PROCID/SLURM_NTASKS, PMI_RANK/PMI_SIZE or
 * OMPI_COMM_WORLD_RANK/OMPI_COMM_WORLD_SIZE); each task works alone and prints
 * its own summary. The exit code is nonzero if any of its files failed.
 *
 * --traj <n> restricts the walk to the top-level directories of src named
 * *.<n>, i.e. the <field>.<traj> directories of one configuration.
 *
 * Usage:
 *   srun -n 256 ./A2AMesonFieldCast --src /path/mf --dst /path/mf_single \
 *       --traj 1540
 */
#include <Hadrons/Global.hpp>
#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

#ifndef HAVE_HDF5
#error "A2AMesonFieldCast requires HDF5 support in Grid/Hadrons"
#endif

using namespace Grid;
using namespace Hadrons;
namespace fs = std::filesystem;

namespace A2AMesonFieldCast
{
struct TaskRank
{
    unsigned long rank, size;
};

struct Options
{
    std::string src, dst;
    long long   traj = -1;
};

struct FileStats
{
    unsigned long long datasets = 0, bytesIn = 0, bytesOut = 0;
};

enum class FileStatus { Converted, Skipped };

static constexpr hsize_t contiguousRows = 256;

void require(bool ok, const std::string &message)
{
    if (!ok) throw std::runtime_error(message);
}

bool exists(const std::string &path)
{
    struct stat st;
    if (::lstat(path.c_str(), &st) == 0) return true;
    require(errno == ENOENT, "Cannot inspect '" + path + "': " + std::strerror(errno));
    return false;
}

unsigned long envNumber(const char *name)
{
    const char *value = std::getenv(name);
    require(value != nullptr, std::string("Launcher variable not set: ") + name);
    char *end = nullptr;
    const unsigned long n = std::strtoul(value, &end, 10);
    require(end != value && *end == '\0', std::string("Invalid launcher variable: ") + name);
    return n;
}

TaskRank taskRank(void)
{
    TaskRank r{0, 1};
    if (std::getenv("SLURM_PROCID"))
        r = {envNumber("SLURM_PROCID"), envNumber("SLURM_NTASKS")};
    else if (std::getenv("PMI_RANK"))
        r = {envNumber("PMI_RANK"), envNumber("PMI_SIZE")};
    else if (std::getenv("OMPI_COMM_WORLD_RANK"))
        r = {envNumber("OMPI_COMM_WORLD_RANK"), envNumber("OMPI_COMM_WORLD_SIZE")};
    require(r.size > 0 && r.rank < r.size, "Invalid task rank/size from the launcher");
    return r;
}

Options parseOptions(int argc, char *argv[])
{
    Options opt;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        require(i + 1 < argc, "Missing value for " + arg);
        const std::string value = argv[++i];
        if (arg == "--src")       opt.src = value;
        else if (arg == "--dst")  opt.dst = value;
        else if (arg == "--traj") opt.traj = std::stoll(value);
        else throw std::runtime_error("Unknown option: " + arg);
    }
    require(!opt.src.empty() && !opt.dst.empty(), "--src and --dst are required");
    return opt;
}

bool isInside(const fs::path &inner, const fs::path &outer)
{
    auto o = outer.begin();
    auto i = inner.begin();
    for (; o != outer.end(); ++o, ++i)
    {
        if (i == inner.end() || *i != *o) return false;
    }
    return true;
}

bool isH5File(const fs::directory_entry &e)
{
    return fs::is_regular_file(e.symlink_status()) && e.path().extension() == ".h5";
}

void walk(const fs::path &dir, const fs::path &root, std::vector<std::string> &files)
{
    for (const auto &e: fs::recursive_directory_iterator(dir))
    {
        if (isH5File(e)) files.push_back(e.path().lexically_relative(root).string());
    }
}

// Relative paths under root, sorted so every task sees the same list.
std::vector<std::string> listFiles(const fs::path &root, long long traj)
{
    std::vector<std::string> files;
    if (traj < 0)
    {
        walk(root, root, files);
    }
    else
    {
        const std::string suffix = "." + std::to_string(traj);
        for (const auto &e: fs::directory_iterator(root))
        {
            const std::string name = e.path().filename().string();
            if (fs::is_directory(e.symlink_status()) && name.size() > suffix.size()
                && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
            {
                walk(e.path(), root, files);
            }
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

// Size of re/im if type is a complex compound, 0 otherwise.
std::size_t complexRealBytes(const H5NS::DataSet &data)
{
    if (data.getDataType().getClass() != H5T_COMPOUND) return 0;
    auto ct = data.getCompType();
    if (ct.getNmembers() != 2) return 0;
    const int re = ct.getMemberIndex("re"), im = ct.getMemberIndex("im");
    if (re < 0 || im < 0) return 0;
    auto rtype = ct.getMemberDataType(static_cast<unsigned int>(re));
    auto itype = ct.getMemberDataType(static_cast<unsigned int>(im));
    const std::size_t bytes = rtype.getSize();
    if (rtype.getClass() != H5T_FLOAT || itype.getClass() != H5T_FLOAT
        || itype.getSize() != bytes || (bytes != 4 && bytes != 8)) return 0;
    return bytes;
}

// From A2AMesonFieldSelectHits: copy attributes of any type, including VL.
void copyAttributes(const H5NS::H5Object &source, H5NS::H5Object &target)
{
    for (int i = 0; i < source.getNumAttrs(); ++i)
    {
        auto attr  = source.openAttribute(static_cast<unsigned int>(i));
        auto type  = attr.getDataType();
        auto space = attr.getSpace();
        require(H5Tdetect_class(type.getId(), H5T_REFERENCE) == 0,
                "HDF5 object/region-reference attributes are not supported");
        auto out = target.createAttribute(attr.getName(), type, space);
        const hssize_t points = space.getSimpleExtentNpoints();
        require(points >= 0, "Invalid attribute dataspace");
        if (points == 0) continue;
        std::vector<unsigned char> buf(static_cast<std::size_t>(points)*type.getSize());
        attr.read(type, buf.data());
        try { out.write(type, buf.data()); }
        catch (...)
        {
            H5Dvlen_reclaim(type.getId(), space.getId(), H5P_DEFAULT, buf.data());
            throw;
        }
        require(H5Dvlen_reclaim(type.getId(), space.getId(), H5P_DEFAULT, buf.data()) >= 0,
                "Could not release variable-length attribute memory");
    }
}

// Rows of whole chunks: (bt, br) along (t, i), all columns.
void bandShape(const H5NS::DSetCreatPropList &plist, const hsize_t dim[3],
               hsize_t &bt, hsize_t &br)
{
    bt = 1;
    br = std::min(contiguousRows, dim[1]);
    if (plist.getLayout() == H5D_CHUNKED)
    {
        hsize_t chunk[3];
        require(plist.getChunk(3, chunk) == 3, "Chunk rank differs from dataset rank");
        bt = chunk[0];
        br = chunk[1];
    }
}

void castDataSet(const H5NS::DataSet &in, H5NS::Group &outGroup, const std::string &name,
                 const std::string &path, FileStats &stats)
{
    auto space = in.getSpace();
    require(space.getSimpleExtentNdims() == 3, "Expected a (t, i, j) complex dataset: " + path);
    hsize_t dim[3];
    space.getSimpleExtentDims(dim);
    auto plist = in.getCreatePlist();
    auto out   = outGroup.createDataSet(name, Hdf5Type<ComplexF>::type(), space, plist);
    copyAttributes(in, out);

    hsize_t bt, br;
    bandShape(plist, dim, bt, br);
    std::vector<ComplexD> d(static_cast<std::size_t>(bt*br*dim[2]));
    std::vector<ComplexF> f(d.size());
    const double largest = static_cast<double>(std::numeric_limits<float>::max());

    for (hsize_t t = 0; t < dim[0]; t += bt)
    for (hsize_t i = 0; i < dim[1]; i += br)
    {
        const hsize_t count[3]  = {std::min(bt, dim[0] - t), std::min(br, dim[1] - i), dim[2]};
        const hsize_t offset[3] = {t, i, 0};
        const std::size_t n     = static_cast<std::size_t>(count[0]*count[1]*count[2]);
        H5NS::DataSpace memSpace(3, count), fileSpace(space);
        fileSpace.selectHyperslab(H5S_SELECT_SET, count, offset);
        in.read(d.data(), Hdf5Type<ComplexD>::type(), memSpace, fileSpace);
        for (std::size_t k = 0; k < n; ++k)
        {
            const double re = d[k].real(), im = d[k].imag();
            require(std::isfinite(re) && std::isfinite(im)
                    && std::abs(re) <= largest && std::abs(im) <= largest,
                    "Non-finite or out-of-float-range value at t=" + std::to_string(t)
                    + " in " + path);
            f[k] = ComplexF(static_cast<float>(re), static_cast<float>(im));
        }
        out.write(f.data(), Hdf5Type<ComplexF>::type(), memSpace, fileSpace);
    }
    const unsigned long long elements = dim[0]*dim[1]*dim[2];
    stats.datasets++;
    stats.bytesIn  += elements*sizeof(ComplexD);
    stats.bytesOut += elements*sizeof(ComplexF);
}

void copyGroup(const H5NS::Group &source, H5NS::Group &target, const std::string &path,
               FileStats &stats, unsigned int depth = 0)
{
    require(depth < 64, "HDF5 hierarchy too deep (possibly cyclic hard links)");
    copyAttributes(source, target);
    for (hsize_t i = 0; i < source.getNumObjs(); ++i)
    {
        const std::string name      = source.getObjnameByIdx(i);
        const std::string childPath = path + "/" + name;
        H5L_info_t linkInfo;
        require(H5Lget_info(source.getId(), name.c_str(), &linkInfo, H5P_DEFAULT) >= 0,
                "Cannot inspect link: " + childPath);
        require(linkInfo.type == H5L_TYPE_HARD, "Soft/external links are not supported: " + childPath);
        const auto objectType = source.getObjTypeByIdx(i);
        if (objectType == H5G_GROUP)
        {
            auto src = source.openGroup(name);
            auto dst = target.createGroup(name);
            copyGroup(src, dst, childPath, stats, depth + 1);
        }
        else if (objectType == H5G_DATASET)
        {
            auto data = source.openDataSet(name);
            const std::size_t bytes = complexRealBytes(data);
            require(bytes != 4, "Dataset is already single precision: " + childPath);
            if (bytes == 8)
            {
                castDataSet(data, target, name, childPath, stats);
            }
            else
            {
                require(H5Tdetect_class(data.getDataType().getId(), H5T_REFERENCE) == 0,
                        "Reference-valued dataset is not supported: " + childPath);
                require(H5Ocopy(source.getId(), name.c_str(), target.getId(), name.c_str(),
                                H5P_DEFAULT, H5P_DEFAULT) >= 0, "Cannot copy dataset: " + childPath);
            }
        }
        else
        {
            require(objectType == H5G_TYPE, "Unsupported HDF5 object: " + childPath);
            require(H5Ocopy(source.getId(), name.c_str(), target.getId(), name.c_str(),
                            H5P_DEFAULT, H5P_DEFAULT) >= 0, "Cannot copy datatype: " + childPath);
        }
    }
}

FileStatus convertFile(const std::string &input, const std::string &output, FileStats &stats)
{
    if (exists(output)) return FileStatus::Skipped;
    const std::string temporary = output + ".partial";
    bool ownTemporary = false;

    try
    {
        require(Hadrons::mkdir(dirname(output)) == 0,
                "Cannot create directory '" + dirname(output) + "': " + std::strerror(errno));
        if (exists(temporary))
        {
            require(::unlink(temporary.c_str()) == 0,
                    "Cannot remove stale " + temporary + ": " + std::strerror(errno));
        }
        {
            H5NS::H5File in(input, H5F_ACC_RDONLY);
            H5NS::H5File out(temporary, H5F_ACC_EXCL);
            ownTemporary = true;
            {
                auto inRoot  = in.openGroup("/");
                auto outRoot = out.openGroup("/");
                copyGroup(inRoot, outRoot, "", stats);
            }
            require(stats.datasets > 0, "No complex double dataset in " + input);
            out.flush(H5F_SCOPE_GLOBAL);
            out.close();
        }
        if (::link(temporary.c_str(), output.c_str()) != 0)
        {
            const int error = errno;
            ownTemporary = false;
            throw std::runtime_error("Cannot publish " + output + ": " + std::strerror(error)
                                     + ". Completed temporary file retained: " + temporary);
        }
        if (::unlink(temporary.c_str()) != 0)
        {
            std::cerr << "Warning: output is complete, but could not remove " << temporary << std::endl;
        }
        ownTemporary = false;
    }
    catch (...)
    {
        if (ownTemporary) ::unlink(temporary.c_str());
        throw;
    }
    return FileStatus::Converted;
}
} // namespace A2AMesonFieldCast

int main(int argc, char *argv[])
{
    using namespace A2AMesonFieldCast;

    H5NS::Exception::dontPrint();
    TaskRank task;
    Options  opt;
    fs::path src, dst;
    std::vector<std::string> files;

    try
    {
        task = taskRank();
        opt  = parseOptions(argc, argv);
        require(fs::is_directory(opt.src), "--src is not a directory: " + opt.src);
        src = fs::canonical(opt.src);
        dst = fs::weakly_canonical(opt.dst);
        require(!isInside(dst, src) && !isInside(src, dst),
                "--src and --dst must be disjoint trees");
        require(Hadrons::mkdir(opt.dst) == 0, "Cannot create --dst: " + opt.dst);
        dst = fs::canonical(opt.dst);
        require(!isInside(dst, src) && !isInside(src, dst),
                "--src and --dst must be disjoint trees");
        files = listFiles(src, opt.traj);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "usage: " << argv[0] << " --src <root> --dst <root> [--traj <n>]" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string tag = "[task " + std::to_string(task.rank) + "/"
                            + std::to_string(task.size) + "] ";
    if (task.rank == 0)
    {
        std::cout << tag << files.size() << " file(s) under " << src << " -> " << dst << std::endl;
    }

    unsigned long nConverted = 0, nSkipped = 0, nFailed = 0;
    unsigned long long bytesIn = 0, bytesOut = 0;
    const auto begin = std::chrono::steady_clock::now();

    for (std::size_t k = task.rank; k < files.size(); k += task.size)
    {
        const std::string input  = (src/files[k]).string();
        const std::string output = (dst/files[k]).string();
        const auto        start  = std::chrono::steady_clock::now();
        FileStats         stats;

        try
        {
            if (convertFile(input, output, stats) == FileStatus::Skipped)
            {
                nSkipped++;
                std::cout << tag << "SKIP " << files[k] << " (exists)" << std::endl;
                continue;
            }
            const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            nConverted++;
            bytesIn  += stats.bytesIn;
            bytesOut += stats.bytesOut;
            std::cout << tag << "OK   " << files[k] << " " << stats.datasets << " dataset(s), "
                      << stats.bytesIn << " -> " << stats.bytesOut << " bytes, " << s << " s"
                      << std::endl;
        }
        catch (const H5NS::Exception &e)
        {
            nFailed++;
            std::cout << tag << "FAIL " << files[k] << " HDF5 error: " << e.getDetailMsg() << std::endl;
        }
        catch (const std::exception &e)
        {
            nFailed++;
            std::cout << tag << "FAIL " << files[k] << " " << e.what() << std::endl;
        }
    }

    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    std::cout << tag << "SUMMARY converted " << nConverted << ", skipped " << nSkipped
              << ", failed " << nFailed << ", payload " << bytesIn << " -> " << bytesOut
              << " bytes, " << seconds << " s" << std::endl;

    return (nFailed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
