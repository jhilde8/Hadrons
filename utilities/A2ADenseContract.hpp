/*
 * A2ADenseContract.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Shared layer of the dense-W meson-field tools: axis descriptions, connection
 * plans, the F/D single-file/timeSliceIO reader and the timeSliceIO writer.
 *
 * Destined for Hadrons/A2AMatrix.hpp once it has settled; it lives under
 * utilities/ until then so that changing it rebuilds no library code.
 */
#ifndef Hadrons_A2ADenseContract_hpp_
#define Hadrons_A2ADenseContract_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <regex>
#include <set>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

#ifndef HAVE_HDF5
#error "A2ADenseContract.hpp requires HDF5 support in Grid/Hadrons"
#endif

BEGIN_HADRONS_NAMESPACE

namespace A2ADense
{

const std::size_t nSc = Ns*Nc;

using Matrix = A2AMatrix<ComplexD>;

/******************************************************************************
 *                              XML description                               *
 ******************************************************************************/
class TrajRange: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                    unsigned int, start,
                                    unsigned int, end,
                                    unsigned int, step);
};

// One axis of a meson field. A dense axis is a compact W leg,
//
//     dense slot    d = nLow + h*nSc + sc            size nLow + nHit*nSc
//
// any other axis is a fully expanded mode index,
//
//     expanded mode e = nLow + (h*nt + t)*nSc + sc   size nLow + nHit*nt*nSc
//
// timeOffset is added to the slice a field is taken at to give the time of
// the W behind a dense axis. It is 0 on every base field; a product's axes
// inherit the offsets of the inputs they come from.
class AxisPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(AxisPar,
                                    bool,         dense,
                                    unsigned int, nLow,
                                    unsigned int, timeOffset);
};

// One meson field, base or product, in every tool. dataset is the HDF5 group
// (the ioname); the matrix is <dataset>/a2aMatrix. With timeSliceIO, file is
// the stem <x>.h5 and slice t is read from <x>.tNNNN.h5.
class MatrixPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MatrixPar,
                                    std::string, name,
                                    std::string, file,
                                    std::string, dataset,
                                    bool,        timeSliceIO,
                                    AxisPar,     row,
                                    AxisPar,     col);
};

// Written next to the matrix of every product slice, and checked against the
// XML description whenever a product is loaded.
class ProductMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ProductMetadata,
                                    AxisPar,      row,
                                    AxisPar,      col,
                                    MatrixPar,    left,
                                    MatrixPar,    right,
                                    unsigned int, deltaL,
                                    unsigned int, deltaR,
                                    unsigned int, nt);
};

const std::string productMetadataName = "A2ADenseProduct";

/******************************************************************************
 *                             Axes and connections                           *
 ******************************************************************************/
struct Axis
{
    bool        dense;
    std::size_t nLow, nHit, timeOffset, size;
};

// nHit is inferred from the size, so a 1-hit AMA field and an 8-hit sloppy
// one need no separate flag.
inline Axis makeAxis(const AxisPar &par, const std::size_t size,
                     const std::size_t nt, const std::string &what)
{
    const std::size_t block = par.dense ? nSc : nSc*nt;

    if ((size < par.nLow) or ((size - par.nLow) % block != 0))
    {
        HADRONS_ERROR(Size, what + ": size " + std::to_string(size)
                      + " is not nLow + nHit*" + std::to_string(block)
                      + " for nLow = " + std::to_string(par.nLow)
                      + (par.dense ? " (dense)" : " (expanded)"));
    }
    if (par.timeOffset >= nt)
    {
        HADRONS_ERROR(Range, what + ": timeOffset "
                      + std::to_string(par.timeOffset) + " >= nt");
    }

    return {par.dense, par.nLow, (size - par.nLow)/block, par.timeOffset, size};
}

// Which operand of a product or trace has its contracted axis gathered: the
// expanded one, whenever the other side is dense.
enum class Gather {none, left, right};

struct Connection
{
    Gather      gather;
    std::size_t nLow, nHit, nt, timeOffset, active;
};

// The edge left.col . right.row. timeOffset is that of the dense axis, whose
// W time selects the expanded side; active is the summed dimension.
inline Connection makeConnection(const Axis &left, const Axis &right,
                                 const std::size_t nt, const std::string &what)
{
    if (left.dense and right.dense)
    {
        HADRONS_ERROR(Argument, what + ": both contracted axes are dense W");
    }
    if ((left.nLow != right.nLow) or (left.nHit != right.nHit))
    {
        HADRONS_ERROR(Size, what + ": contracted axes differ (nLow "
                      + std::to_string(left.nLow) + "/" + std::to_string(right.nLow)
                      + ", nHit " + std::to_string(left.nHit) + "/"
                      + std::to_string(right.nHit) + ")");
    }
    if (!left.dense and !right.dense)
    {
        return {Gather::none, left.nLow, left.nHit, nt, 0, left.size};
    }

    const Axis &dense = left.dense ? left : right;

    return {left.dense ? Gather::right : Gather::left,
            dense.nLow, dense.nHit, nt, dense.timeOffset, dense.size};
}

// Expanded partner of dense slot d when the dense W sits at time t.
inline std::size_t expandIndex(const Connection &c, const std::size_t d,
                               const std::size_t t)
{
    if (d < c.nLow)
    {
        return d;
    }

    const std::size_t h = (d - c.nLow)/nSc, sc = (d - c.nLow)%nSc;

    return c.nLow + (h*c.nt + t)*nSc + sc;
}

// Expanded rows (cols) of in that meet a dense axis whose W sits at time t:
// the low block plus one strip of nSc per hit. Out may be row- or
// column-major; the copy runs over the rows of in, which are contiguous.
template <typename Out, typename In>
inline void gatherRows(Out &out, const In &in, const Connection &c,
                       const std::size_t t)
{
    out.resize(c.active, in.cols());
    thread_for(r, c.active,
    {
        out.row(r) = in.row(expandIndex(c, r, t));
    });
}

template <typename Out, typename In>
inline void gatherCols(Out &out, const In &in, const Connection &c,
                       const std::size_t t)
{
    const Eigen::Index sc = nSc, low = c.nLow, nr = in.rows();

    out.resize(in.rows(), c.active);
    thread_for(r, nr,
    {
        out.row(r).head(low) = in.row(r).head(low);
        for (std::size_t h = 0; h < c.nHit; ++h)
        {
            out.row(r).segment(low + h*sc, sc) =
                in.row(r).segment(low + (h*c.nt + t)*sc, sc);
        }
    });
}

// out = a*b over the connection c, with a and b taken at slices ta and tb.
// The gathered operand goes to buf, which is reused across calls.
inline void multiply(Matrix &out, const Matrix &a, const std::size_t ta,
                     const Matrix &b, const std::size_t tb, const Connection &c,
                     Matrix &buf)
{
    const Matrix *ap = &a, *bp = &b;

    if (c.gather == Gather::left)
    {
        gatherCols(buf, a, c, (tb + c.timeOffset) % c.nt);
        ap = &buf;
    }
    else if (c.gather == Gather::right)
    {
        gatherRows(buf, b, c, (ta + c.timeOffset) % c.nt);
        bp = &buf;
    }
    out.noalias() = (*ap)*(*bp);
}

using MatrixTr = A2AMatrixTr<ComplexD>;

struct TraceBuffers
{
    Matrix   a, tmp;
    MatrixTr b;
};

// out = in with its axes gathered where their edges say so. An edge is always
// left.col . right.row, so a matrix's rows are the right side of rowEdge and
// its cols the left side of colEdge. tRow/tCol are the partner slices on the
// other end of each edge; both gathered (VV against WW) goes through tmp.
template <typename Out>
inline void select(Out &out, const Matrix &in,
                   const Connection &rowEdge, const std::size_t tRow,
                   const Connection &colEdge, const std::size_t tCol, Matrix &tmp)
{
    const bool        rows = (rowEdge.gather == Gather::right);
    const bool        cols = (colEdge.gather == Gather::left);
    const std::size_t tr   = (tRow + rowEdge.timeOffset) % rowEdge.nt;
    const std::size_t tc   = (tCol + colEdge.timeOffset) % colEdge.nt;

    if (rows and cols)
    {
        gatherRows(tmp, in, rowEdge, tr);
        gatherCols(out, tmp, colEdge, tc);
    }
    else if (rows)
    {
        gatherRows(out, in, rowEdge, tr);
    }
    else if (cols)
    {
        gatherCols(out, in, colEdge, tc);
    }
    else
    {
        out = in;
    }
}

// in reduced as select does, before it enters any product: an expanded axis
// that meets one W time is gathered first, so no product carries columns a
// later gather would drop. Returns in itself, uncopied, when nothing is
// gathered. An open axis (one closed later at many times) takes openEdge.
inline const Matrix &gathered(Matrix &out, const Matrix &in,
                              const Connection &rowEdge, const std::size_t tRow,
                              const Connection &colEdge, const std::size_t tCol,
                              Matrix &tmp)
{
    if ((rowEdge.gather != Gather::right) and (colEdge.gather != Gather::left))
    {
        return in;
    }
    select(out, in, rowEdge, tRow, colEdge, tCol, tmp);

    return out;
}

const Connection openEdge{Gather::none, 0, 0, 1, 0, 0};

// acc += tr(a*b) over inner = a.col . b.row and closing = b.col . a.row.
// aRowTime and aColTime are the slices behind a's row and col axes, which
// differ when a is a running product; b is one field taken at tb.
inline void accTrace(ComplexD &acc, const Matrix &a, const std::size_t aRowTime,
                     const std::size_t aColTime, const Matrix &b, const std::size_t tb,
                     const Connection &inner, const Connection &closing,
                     TraceBuffers &buf)
{
    const Matrix *ap = &a;

    if ((closing.gather == Gather::right) or (inner.gather == Gather::left))
    {
        select(buf.a, a, closing, tb, inner, tb, buf.tmp);
        ap = &buf.a;
    }
    // b always goes to the column-major buffer, a plain copy if nothing is gathered.
    select(buf.b, b, inner, aColTime, closing, aRowTime, buf.tmp);
    A2AContraction::accTrMul(acc, *ap, buf.b);
}

// acc += tr(m) for m taken at slice t, closed on itself through edge =
// m.col . m.row. The sum is O(rows), so it runs serially.
inline void accSelfTrace(ComplexD &acc, const Matrix &m, const std::size_t t,
                         const Connection &edge)
{
    const std::size_t tw = (t + edge.timeOffset) % edge.nt;

    if (edge.gather == Gather::left)
    {
        for (std::size_t d = 0; d < edge.active; ++d)
        {
            acc += m(d, expandIndex(edge, d, tw));
        }
    }
    else if (edge.gather == Gather::right)
    {
        for (std::size_t d = 0; d < edge.active; ++d)
        {
            acc += m(expandIndex(edge, d, tw), d);
        }
    }
    else
    {
        acc += m.diagonal().sum();
    }
}

inline double seconds(const std::chrono::steady_clock::time_point &start)
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// "0..3 7" -> {0, 1, 2, 3, 7}; every value must be < nt.
inline std::set<unsigned int> parseTimeRange(const std::string &str, const unsigned int nt)
{
    const std::regex       rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
    std::set<unsigned int> tSet;

    for (auto &s: strToVec<std::string>(str))
    {
        std::smatch sm;

        if (!std::regex_match(s, sm, rex))
        {
            HADRONS_ERROR(Argument, "invalid time expression '" + s + "'");
        }

        const unsigned long a = std::stoul(sm[1].matched ? sm[1].str() : sm[3].str());
        const unsigned long b = sm[1].matched ? a : std::stoul(sm[4].str());

        if ((a > b) or (b >= nt))
        {
            HADRONS_ERROR(Range, "invalid time range '" + s + "' for nt = "
                          + std::to_string(nt));
        }
        for (unsigned long t = a; t <= b; ++t)
        {
            tSet.insert(t);
        }
    }

    return tSet;
}

/******************************************************************************
 *                                    I/O                                     *
 ******************************************************************************/
inline std::string sliceFilename(const std::string &stem, const unsigned int t)
{
    const std::string  ext = ".h5";
    std::ostringstream ss;

    if ((stem.size() < ext.size())
        or (stem.compare(stem.size() - ext.size(), ext.size(), ext) != 0))
    {
        HADRONS_ERROR(Io, "timeSliceIO file '" + stem + "' does not end in '.h5'");
    }
    ss << stem.substr(0, stem.size() - ext.size()) << ".t"
       << std::setfill('0') << std::setw(4) << t << ext;

    return ss.str();
}

inline bool fileExists(const std::string &path)
{
    struct stat s;

    return (stat(path.c_str(), &s) == 0);
}

// A resident field: every slice held in memory as ComplexD.
struct Field
{
    MatrixPar           par;
    Axis                row, col;
    std::vector<Matrix> slice;
    bool                isProduct{false};
};

// The memory type is always ComplexD, so HDF5 converts a ComplexF file on
// read; the file's own type only has to be a float re/im compound.
inline std::size_t checkComplexType(const H5NS::DataSet &data, const std::string &what)
{
    if (data.getTypeClass() != H5T_COMPOUND)
    {
        HADRONS_ERROR(Io, what + ": matrix is not a re/im compound");
    }

    const auto type = data.getCompType();
    const int  re   = H5Tget_member_index(type.getId(), "re");
    const int  im   = H5Tget_member_index(type.getId(), "im");

    if ((type.getNmembers() != 2) or (re < 0) or (im < 0))
    {
        HADRONS_ERROR(Io, what + ": matrix is not a re/im compound");
    }

    const auto rt = type.getMemberDataType(re), it = type.getMemberDataType(im);

    if ((rt.getClass() != H5T_FLOAT) or (it.getClass() != H5T_FLOAT)
        or (rt.getSize() != it.getSize())
        or ((rt.getSize() != 4) and (rt.getSize() != 8)))
    {
        HADRONS_ERROR(Io, what + ": re/im are not both float or both double");
    }

    return rt.getSize();
}

struct SliceDims
{
    hsize_t nt, ni, nj;
};

inline SliceDims readDims(const H5NS::DataSet &data, const std::string &what)
{
    hsize_t d[3];
    auto    space = data.getSpace();

    if (space.getSimpleExtentNdims() != 3)
    {
        HADRONS_ERROR(Io, what + ": matrix is not (t, row, col)");
    }
    space.getSimpleExtentDims(d);

    return {d[0], d[1], d[2]};
}

inline void readSlice(Matrix &m, const H5NS::DataSet &data, const hsize_t t,
                      const hsize_t ni, const hsize_t nj, const std::string &what)
{
    const hsize_t   count[3] = {1, ni, nj}, offset[3] = {t, 0, 0};
    const hsize_t   memCount[2] = {ni, nj};
    H5NS::DataSpace space = data.getSpace(), mem(2, memCount);

    m.resize(ni, nj);
    space.selectHyperslab(H5S_SELECT_SET, count, offset);
    data.read(m.data(), Hdf5Type<ComplexD>::type(), mem, space);
    if (!m.allFinite())
    {
        HADRONS_ERROR(Io, what + ": non-finite entry at t = " + std::to_string(t));
    }
}

// A product must carry metadata matching its XML axes; a base field carries
// none and is checked on its dimensions alone.
inline bool checkProductMetadata(const std::string &filename, const MatrixPar &par)
{
    bool hasMetadata;

    {
        H5NS::H5File file(filename, H5F_ACC_RDONLY);
        auto         group = file.openGroup(par.dataset);

        hasMetadata = (H5Lexists(group.getId(), productMetadataName.c_str(),
                                 H5P_DEFAULT) > 0);
    }
    if (hasMetadata)
    {
        Hdf5Reader      reader(filename);
        ProductMetadata md;

        push(reader, par.dataset);
        read(reader, productMetadataName, md);
        pop(reader);
        if (!(md.row == par.row) or !(md.col == par.col))
        {
            HADRONS_ERROR(Definition, "'" + par.name + "': XML axes differ from the "
                          "product metadata in '" + filename + "'");
        }
    }

    return hasMetadata;
}

// One field opened for slice-by-slice reading: shape, type and product
// metadata are checked once by openSliceReader. Each read reopens its file
// (the slice file with timeSliceIO), so no handle is held between reads.
struct SliceReader
{
    MatrixPar    par;
    unsigned int nt;
    Axis         row, col;
    bool         isProduct;

    void read(Matrix &m, const unsigned int t) const
    {
        const std::string fn   = par.timeSliceIO ? sliceFilename(par.file, t) : par.file;
        H5NS::H5File      file(fn, H5F_ACC_RDONLY);
        auto              data = file.openGroup(par.dataset).openDataSet(HADRONS_A2AM_NAME);
        const SliceDims   d    = readDims(data, fn);
        const hsize_t     nts  = par.timeSliceIO ? 1 : nt;

        if ((d.nt != nts) or (d.ni != row.size) or (d.nj != col.size))
        {
            HADRONS_ERROR(Size, fn + ": matrix is " + std::to_string(d.nt) + "x"
                          + std::to_string(d.ni) + "x" + std::to_string(d.nj)
                          + ", expected " + std::to_string(nts) + "x"
                          + std::to_string(row.size) + "x" + std::to_string(col.size));
        }
        checkComplexType(data, fn);
        readSlice(m, data, par.timeSliceIO ? 0 : t, d.ni, d.nj, fn);
    }
};

// par.file must already be trajectory substituted.
inline SliceReader openSliceReader(const MatrixPar &par, const unsigned int nt)
{
    const std::string fn = par.timeSliceIO ? sliceFilename(par.file, 0) : par.file;
    const hsize_t     nts = par.timeSliceIO ? 1 : nt;
    SliceReader       s;
    SliceDims         d;
    std::size_t       realBytes;

    {
        H5NS::H5File file(fn, H5F_ACC_RDONLY);
        auto         data = file.openGroup(par.dataset).openDataSet(HADRONS_A2AM_NAME);

        d = readDims(data, fn);
        if (d.nt != nts)
        {
            HADRONS_ERROR(Size, fn + ": " + std::to_string(d.nt) + " timeslices, expected "
                          + std::to_string(nts));
        }
        realBytes = checkComplexType(data, fn);
    }
    s.par       = par;
    s.nt        = nt;
    s.isProduct = checkProductMetadata(fn, par);
    s.row       = makeAxis(par.row, d.ni, nt, "'" + par.name + "' rows");
    s.col       = makeAxis(par.col, d.nj, nt, "'" + par.name + "' cols");
    std::cout << "  '" << par.name << "' " << nt << "x" << d.ni << "x" << d.nj
              << " Complex" << ((realBytes == 4) ? "F" : "D")
              << (s.isProduct ? " (product)" : "") << std::endl;

    return s;
}

// Loads every slice of one field.
inline Field loadField(const MatrixPar &par, const unsigned int nt)
{
    const SliceReader reader = openSliceReader(par, nt);
    Field             f;

    f.par       = par;
    f.row       = reader.row;
    f.col       = reader.col;
    f.isProduct = reader.isProduct;
    f.slice.resize(nt);
    for (unsigned int t = 0; t < nt; ++t)
    {
        reader.read(f.slice[t], t);
    }

    return f;
}

// Writes one ComplexD slice as a (1, rows, cols) file through a temporary
// renamed into place, so a slice file that exists is complete. Never
// overwrites.
inline void writeSlice(const std::string &filename, const std::string &dataset,
                       const Matrix &m, const ProductMetadata &md)
{
    const hsize_t     chunkSize   = 256;
    const hsize_t     ni = m.rows(), nj = m.cols();
    const hsize_t     dim[3]      = {1, ni, nj};
    const hsize_t     chunk[3]    = {1, std::min(chunkSize, ni), std::min(chunkSize, nj)};
    const hsize_t     memCount[2] = {ni, nj};
    const std::string tmp = filename + ".partial." + std::to_string(getpid());

    if (fileExists(filename))
    {
        HADRONS_ERROR(Io, "output '" + filename + "' already exists");
    }
    {
        Hdf5Writer writer(tmp);

        push(writer, dataset);
        write(writer, productMetadataName, md);

        H5NS::DSetCreatPropList plist;
        H5NS::DataSpace         space(3, dim), mem(2, memCount);

        plist.setChunk(3, chunk);
        plist.setFletcher32();

        auto data = writer.getGroup().createDataSet(HADRONS_A2AM_NAME,
                                                    Hdf5Type<ComplexD>::type(),
                                                    space, plist);

        data.write(m.data(), Hdf5Type<ComplexD>::type(), mem, space);
    }
    if (std::rename(tmp.c_str(), filename.c_str()) != 0)
    {
        HADRONS_ERROR(Io, "cannot rename '" + tmp + "' to '" + filename + "': "
                      + std::strerror(errno));
    }
}

}

END_HADRONS_NAMESPACE

#endif
