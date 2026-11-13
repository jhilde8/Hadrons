/*
 * Test_mf_cast_compare.cpp
 *
 * Checks a single-precision meson-field tree written by A2AMesonFieldCast
 * against the double-precision tree it was converted from. Every
 * <double>/<rel>.h5 must have a <single>/<rel>.h5 with the same HDF5 tree:
 *   - the same groups, datasets, named datatypes and attributes, with equal
 *     types, extents and values (variable-length strings compared as
 *     strings, everything else byte for byte);
 *   - every complex double dataset stored as Hdf5Type<ComplexF> with the same
 *     dims, layout, chunk dims and filters, and every element equal to the
 *     round-to-nearest float cast of the double one (re and im separately).
 *     This is exact, not a tolerance. The worst relative difference over
 *     components in the normal float range (|d| >= FLT_MIN) is printed as
 *     well, and is at most 2^-24; nonzero components below that range lose
 *     relative precision (down to 0) and are counted separately.
 * Task 0 also reports *.h5 files under <single> with no double counterpart.
 *
 * As in A2AMesonFieldCast, the file list is split round-robin over the
 * launcher's tasks without MPI, each task prints its own summary, and the
 * exit code is nonzero if any of its files failed. --traj restricts both
 * walks to the top-level directories named *.<traj>.
 *
 * Usage:
 *   srun -n 256 ./Test_mf_cast_compare --double /path/mf --single /path/mf_single \
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
#include <set>
#include <stdexcept>

#ifndef HAVE_HDF5
#error "Test_mf_cast_compare requires HDF5 support in Grid/Hadrons"
#endif

using namespace Grid;
using namespace Hadrons;
namespace fs = std::filesystem;

struct TaskRank
{
    unsigned long rank, size;
};

struct Options
{
    std::string doubleRoot, singleRoot;
    long long   traj = -1;
};

struct CompareStats
{
    unsigned long long datasets = 0, elements = 0, mismatches = 0, belowNormal = 0;
    double             maxRel   = 0.;
    std::string        firstMismatch;
};

static constexpr hsize_t contiguousRows = 256;

static void require(bool ok, const std::string &message)
{
    if (!ok) throw std::runtime_error(message);
}

static unsigned long envNumber(const char *name)
{
    const char *value = std::getenv(name);
    require(value != nullptr, std::string("Launcher variable not set: ") + name);
    char *end = nullptr;
    const unsigned long n = std::strtoul(value, &end, 10);
    require(end != value && *end == '\0', std::string("Invalid launcher variable: ") + name);
    return n;
}

static TaskRank taskRank(void)
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

static Options parseOptions(int argc, char *argv[])
{
    Options opt;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        require(i + 1 < argc, "Missing value for " + arg);
        const std::string value = argv[++i];
        if (arg == "--double")      opt.doubleRoot = value;
        else if (arg == "--single") opt.singleRoot = value;
        else if (arg == "--traj")   opt.traj = std::stoll(value);
        else throw std::runtime_error("Unknown option: " + arg);
    }
    require(!opt.doubleRoot.empty() && !opt.singleRoot.empty(), "--double and --single are required");
    return opt;
}

static void walk(const fs::path &dir, const fs::path &root, std::vector<std::string> &files)
{
    for (const auto &e: fs::recursive_directory_iterator(dir))
    {
        if (fs::is_regular_file(e.symlink_status()) && e.path().extension() == ".h5")
        {
            files.push_back(e.path().lexically_relative(root).string());
        }
    }
}

static std::vector<std::string> listFiles(const fs::path &root, long long traj)
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

static std::size_t complexRealBytes(const H5NS::DataSet &data)
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

static bool typesEqual(const H5NS::DataType &a, const H5NS::DataType &b)
{
    return H5Tequal(a.getId(), b.getId()) > 0;
}

static void requireSameExtent(const H5NS::DataSpace &a, const H5NS::DataSpace &b,
                              const std::string &path)
{
    require(H5Sextent_equal(a.getId(), b.getId()) > 0, "Extent differs: " + path);
}

static void readValue(const H5NS::Attribute &a, const H5NS::DataType &t, void *buf) { a.read(t, buf); }
static void readValue(const H5NS::DataSet &d, const H5NS::DataType &t, void *buf) { d.read(buf, t); }
static void readValue(const H5NS::Attribute &a, const H5NS::DataType &t, std::string &s) { a.read(t, s); }
static void readValue(const H5NS::DataSet &d, const H5NS::DataType &t, std::string &s) { d.read(s, t); }

// Raw values of a fixed-size type; variable-length strings as strings.
template <typename Obj>
static void requireSameValues(const Obj &a, const Obj &b, const H5NS::DataType &type,
                              hssize_t points, const std::string &path)
{
    require(points >= 0, "Invalid dataspace: " + path);
    if (H5Tis_variable_str(type.getId()) > 0)
    {
        require(points == 1, "Non-scalar variable-length string is not supported: " + path);
        std::string sa, sb;
        readValue(a, type, sa);
        readValue(b, type, sb);
        require(sa == sb, "String value differs: " + path);
        return;
    }
    require(H5Tdetect_class(type.getId(), H5T_VLEN) == 0
            && H5Tdetect_class(type.getId(), H5T_REFERENCE) == 0,
            "Variable-length or reference values are not supported: " + path);
    const std::size_t bytes = static_cast<std::size_t>(points)*type.getSize();
    std::vector<unsigned char> va(bytes), vb(bytes);
    if (bytes == 0) return;
    readValue(a, type, va.data());
    readValue(b, type, vb.data());
    require(std::memcmp(va.data(), vb.data(), bytes) == 0, "Value differs: " + path);
}

static void compareAttributes(const H5NS::H5Object &a, const H5NS::H5Object &b,
                              const std::string &path)
{
    require(a.getNumAttrs() == b.getNumAttrs(), "Attribute count differs: " + path);
    for (int i = 0; i < a.getNumAttrs(); ++i)
    {
        auto aa = a.openAttribute(static_cast<unsigned int>(i));
        const std::string name = aa.getName(), attrPath = path + "@" + name;
        require(b.attrExists(name), "Attribute missing: " + attrPath);
        auto ab = b.openAttribute(name);
        auto type = aa.getDataType();
        require(typesEqual(type, ab.getDataType()), "Attribute type differs: " + attrPath);
        requireSameExtent(aa.getSpace(), ab.getSpace(), attrPath);
        requireSameValues(aa, ab, type, aa.getSpace().getSimpleExtentNpoints(), attrPath);
    }
}

static void comparePlists(const H5NS::DSetCreatPropList &a, const H5NS::DSetCreatPropList &b,
                          const std::string &path)
{
    require(a.getLayout() == b.getLayout(), "Layout differs: " + path);
    if (a.getLayout() == H5D_CHUNKED)
    {
        hsize_t ca[3], cb[3];
        require(a.getChunk(3, ca) == 3 && b.getChunk(3, cb) == 3, "Chunk rank is not 3: " + path);
        require(std::equal(ca, ca + 3, cb), "Chunk dims differ: " + path);
    }
    const int n = a.getNfilters();
    require(n == b.getNfilters(), "Filter count differs: " + path);
    for (int k = 0; k < n; ++k)
    {
        unsigned int flags, cd[16], config;
        std::size_t  nd = 16;
        char         name[64];
        const H5Z_filter_t fa = H5Pget_filter2(a.getId(), k, &flags, &nd, cd, sizeof(name), name, &config);
        nd = 16;
        const H5Z_filter_t fb = H5Pget_filter2(b.getId(), k, &flags, &nd, cd, sizeof(name), name, &config);
        require(fa >= 0 && fa == fb, "Filter " + std::to_string(k) + " differs: " + path);
    }
}

static void compareComplex(const H5NS::DataSet &d, const H5NS::DataSet &s,
                           const std::string &path, CompareStats &stats)
{
    require(typesEqual(s.getDataType(), Hdf5Type<ComplexF>::type()),
            "Single dataset is not Hdf5Type<ComplexF>: " + path);
    auto space = d.getSpace();
    require(space.getSimpleExtentNdims() == 3, "Expected a (t, i, j) complex dataset: " + path);
    hsize_t dim[3];
    space.getSimpleExtentDims(dim);
    auto plist = d.getCreatePlist();
    comparePlists(plist, s.getCreatePlist(), path);

    hsize_t bt = 1, br = std::min(contiguousRows, dim[1]);
    if (plist.getLayout() == H5D_CHUNKED)
    {
        hsize_t chunk[3];
        plist.getChunk(3, chunk);
        bt = chunk[0];
        br = chunk[1];
    }
    std::vector<ComplexD> vd(static_cast<std::size_t>(bt*br*dim[2]));
    std::vector<ComplexF> vs(vd.size());

    for (hsize_t t = 0; t < dim[0]; t += bt)
    for (hsize_t i = 0; i < dim[1]; i += br)
    {
        const hsize_t count[3]  = {std::min(bt, dim[0] - t), std::min(br, dim[1] - i), dim[2]};
        const hsize_t offset[3] = {t, i, 0};
        const std::size_t n     = static_cast<std::size_t>(count[0]*count[1]*count[2]);
        H5NS::DataSpace memSpace(3, count), fileSpace(space);
        fileSpace.selectHyperslab(H5S_SELECT_SET, count, offset);
        d.read(vd.data(), Hdf5Type<ComplexD>::type(), memSpace, fileSpace);
        s.read(vs.data(), Hdf5Type<ComplexF>::type(), memSpace, fileSpace);
        for (std::size_t k = 0; k < n; ++k)
        {
            const double dc[2] = {vd[k].real(), vd[k].imag()};
            const float  sc[2] = {vs[k].real(), vs[k].imag()};
            for (int c = 0; c < 2; ++c)
            {
                if (sc[c] != static_cast<float>(dc[c]))
                {
                    if (stats.mismatches == 0)
                    {
                        const hsize_t r = k/dim[2], j = k%dim[2];
                        stats.firstMismatch = path + " at (" + std::to_string(t + r/count[1]) + ", "
                                              + std::to_string(i + r%count[1]) + ", "
                                              + std::to_string(j) + ")";
                    }
                    stats.mismatches++;
                }
                if (std::abs(dc[c]) >= static_cast<double>(std::numeric_limits<float>::min()))
                {
                    stats.maxRel = std::max(stats.maxRel, std::abs(sc[c] - dc[c])/std::abs(dc[c]));
                }
                else if (dc[c] != 0.)
                {
                    stats.belowNormal++;
                }
            }
        }
    }
    stats.datasets++;
    stats.elements += dim[0]*dim[1]*dim[2];
}

static void compareDataSet(const H5NS::DataSet &d, const H5NS::DataSet &s,
                           const std::string &path, CompareStats &stats)
{
    requireSameExtent(d.getSpace(), s.getSpace(), path);
    compareAttributes(d, s, path);
    const std::size_t bytes = complexRealBytes(d);
    require(bytes != 4, "Double-tree dataset is single precision: " + path);
    if (bytes == 8)
    {
        compareComplex(d, s, path, stats);
    }
    else
    {
        auto type = d.getDataType();
        require(typesEqual(type, s.getDataType()), "Dataset type differs: " + path);
        requireSameValues(d, s, type, d.getSpace().getSimpleExtentNpoints(), path);
    }
}

static void compareGroup(const H5NS::Group &d, const H5NS::Group &s, const std::string &path,
                         CompareStats &stats, unsigned int depth = 0)
{
    require(depth < 64, "HDF5 hierarchy too deep (possibly cyclic hard links)");
    compareAttributes(d, s, path.empty() ? "/" : path);
    require(d.getNumObjs() == s.getNumObjs(), "Object count differs: " + path + "/");
    for (hsize_t i = 0; i < d.getNumObjs(); ++i)
    {
        const std::string name = d.getObjnameByIdx(i), childPath = path + "/" + name;
        require(H5Lexists(s.getId(), name.c_str(), H5P_DEFAULT) > 0, "Object missing: " + childPath);
        const H5O_type_t type = d.childObjType(name);
        require(type == s.childObjType(name), "Object kind differs: " + childPath);
        if (type == H5O_TYPE_GROUP)
        {
            compareGroup(d.openGroup(name), s.openGroup(name), childPath, stats, depth + 1);
        }
        else if (type == H5O_TYPE_DATASET)
        {
            compareDataSet(d.openDataSet(name), s.openDataSet(name), childPath, stats);
        }
        else
        {
            require(type == H5O_TYPE_NAMED_DATATYPE, "Unsupported HDF5 object: " + childPath);
            require(typesEqual(d.openDataType(name), s.openDataType(name)),
                    "Named datatype differs: " + childPath);
        }
    }
}

int main(int argc, char *argv[])
{
    H5NS::Exception::dontPrint();
    TaskRank task;
    Options  opt;
    fs::path doubleRoot, singleRoot;
    std::vector<std::string> files;

    try
    {
        task = taskRank();
        opt  = parseOptions(argc, argv);
        require(fs::is_directory(opt.doubleRoot), "--double is not a directory: " + opt.doubleRoot);
        require(fs::is_directory(opt.singleRoot), "--single is not a directory: " + opt.singleRoot);
        doubleRoot = fs::canonical(opt.doubleRoot);
        singleRoot = fs::canonical(opt.singleRoot);
        files      = listFiles(doubleRoot, opt.traj);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "usage: " << argv[0] << " --double <root> --single <root> [--traj <n>]" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string tag = "[task " + std::to_string(task.rank) + "/"
                            + std::to_string(task.size) + "] ";
    unsigned long      nPass = 0, nFail = 0;
    unsigned long long belowNormal = 0;
    double             maxRel = 0.;
    const auto    begin  = std::chrono::steady_clock::now();

    if (task.rank == 0)
    {
        std::cout << tag << files.size() << " file(s) under " << doubleRoot << " vs " << singleRoot
                  << std::endl;
        try
        {
            const std::set<std::string> known(files.begin(), files.end());
            for (const auto &f: listFiles(singleRoot, opt.traj))
            {
                if (!known.count(f))
                {
                    nFail++;
                    std::cout << tag << "FAIL " << f << " has no double counterpart" << std::endl;
                }
            }
        }
        catch (const std::exception &e)
        {
            nFail++;
            std::cout << tag << "FAIL cannot list " << singleRoot << ": " << e.what() << std::endl;
        }
    }

    for (std::size_t k = task.rank; k < files.size(); k += task.size)
    {
        const std::string dFile = (doubleRoot/files[k]).string();
        const std::string sFile = (singleRoot/files[k]).string();
        CompareStats      stats;

        try
        {
            require(fs::is_regular_file(fs::symlink_status(sFile)), "Single file missing");
            H5NS::H5File d(dFile, H5F_ACC_RDONLY), s(sFile, H5F_ACC_RDONLY);
            compareGroup(d.openGroup("/"), s.openGroup("/"), "", stats);
            require(stats.datasets > 0, "No complex double dataset");
            require(stats.mismatches == 0, std::to_string(stats.mismatches)
                    + " component(s) differ from the float cast, first " + stats.firstMismatch);
            nPass++;
            maxRel = std::max(maxRel, stats.maxRel);
            belowNormal += stats.belowNormal;
            std::cout << tag << "PASS " << files[k] << " " << stats.datasets << " dataset(s), "
                      << stats.elements << " elements, max rel diff " << stats.maxRel
                      << ", components below float normal range " << stats.belowNormal << std::endl;
        }
        catch (const H5NS::Exception &e)
        {
            nFail++;
            std::cout << tag << "FAIL " << files[k] << " HDF5 error: " << e.getDetailMsg() << std::endl;
        }
        catch (const std::exception &e)
        {
            nFail++;
            std::cout << tag << "FAIL " << files[k] << " " << e.what() << std::endl;
        }
    }

    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    std::cout << tag << "SUMMARY pass " << nPass << ", fail " << nFail << ", max rel diff "
              << maxRel << " (2^-24 = " << std::ldexp(1., -24) << "), components below float normal range "
              << belowNormal << ", " << seconds << " s"
              << std::endl;

    return (nFail == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
