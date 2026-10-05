/*
 * A2AVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef A2A_Vectors_hpp_
#define A2A_Vectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  High-mode A2A vector construction. Entirely independent of the Schur      *
 *  convention: V is whatever the injected solver returns and W is the raw    *
 *  noise, so nothing here builds or touches a Schur operator. Kept separate  *
 *  from the low-mode classes below so that a high-mode-only module (e.g.     *
 *  MSolver::A2AHighModeVBinned) does not allocate their four red-black 5D    *
 *  scratch fields, which only the low-mode reconstruction needs.             *
 ******************************************************************************/
template <typename FImpl>
class A2AHighModes
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2AHighModes(FMat &action, Solver &solver);
    void makeHighModeV(FermionField &vout, const FermionField &noise);
    void makeHighModeV5D(FermionField &vout_4d, FermionField &vout_5d,
                         const FermionField &noise_5d);
    //////////////////////////////////////////////////////////////////////////
    // Batched forms: one solver call for the whole batch. A solver whose
    // OperatorFunction implements a vector overload then sees every source at
    // once, and the vector SchurRedBlackBase::operator() invokes its guesser
    // once per batch rather than once per source. A solver without that
    // overload still solves the batch one source at a time through
    // OperatorFunction's base-class loop, so these are correct either way.
    // vout_5d is caller-owned scratch, as in the single-source form; the
    // imported 5D sources are local to the call.
    //////////////////////////////////////////////////////////////////////////
    void makeHighModeV(std::vector<FermionField> &vout,
                       const std::vector<FermionField> &noise);
    void makeHighModeV5D(std::vector<FermionField> &vout_4d,
                         std::vector<FermionField> &vout_5d,
                         const std::vector<FermionField> &noise);
    void makeHighModeW(FermionField &wout, const FermionField &noise);
    void makeHighModeW5D(FermionField &vout_5d, FermionField &wout_5d,
                         const FermionField &noise_5d);
protected:
    FMat         &action_;
    Solver       &solver_;
    GridBase     *fGrid_;
    FermionField tmp5_;
};

/******************************************************************************
 *  Shared state for the Schur-convention-specific low-mode constructions     *
 *  below. The convention is a compile-time property of the build             *
 *  (HADRONS_DEFAULT_SCHUR in Global.hpp, one build per ensemble), so this is *
 *  a plain base holding the fields both conventions need -- there is no      *
 *  runtime dispatch and therefore nothing virtual. Modules reach a concrete  *
 *  class through HADRONS_DEFAULT_SCHUR_A2A and never name a convention.      *
 ******************************************************************************/
template <typename FImpl>
class A2ALowModesSchurBase
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2ALowModesSchurBase(FMat &action);
protected:
    FMat         &action_;
    GridBase     *fGrid_, *frbGrid_;
    FermionField tmp5_;
    FermionField src_o_, sol_e_, sol_o_, tmp_;
};

/******************************************************************************
 *              Low-mode V & W, DiagTwo Schur convention                      *
 ******************************************************************************/
template <typename FImpl>
class A2ALowModesSchurDiagTwo : public A2ALowModesSchurBase<FImpl>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2ALowModesSchurDiagTwo(FMat &action);
    void makeLowModeV(FermionField &vout,
                      const FermionField &evec, const Real &eval);
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d,
                        const FermionField &evec, const Real &eval);
    void makeLowModeW(FermionField &wout,
                      const FermionField &evec, const Real &eval);
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d,
                        const FermionField &evec, const Real &eval);
    // Exposes the Schur operator used internally by makeLowModeW/op_, so
    // callers (e.g. A2ALowModeCoarseBinned's in-program eigenvector check) can
    // validate evec_i against the exact same operator this class uses to build
    // V/W rather than building a second one themselves.
    SchurOperatorBase<FermionField>& op(void);
private:
    using A2ALowModesSchurBase<FImpl>::action_;
    using A2ALowModesSchurBase<FImpl>::fGrid_;
    using A2ALowModesSchurBase<FImpl>::frbGrid_;
    using A2ALowModesSchurBase<FImpl>::tmp5_;
    using A2ALowModesSchurBase<FImpl>::src_o_;
    using A2ALowModesSchurBase<FImpl>::sol_e_;
    using A2ALowModesSchurBase<FImpl>::sol_o_;
    using A2ALowModesSchurBase<FImpl>::tmp_;
    SchurDiagTwoOperator<FMat, FermionField> op_;
};

template <typename FImpl>
class A2ALowModesSchurDiagOne : public A2ALowModesSchurBase<FImpl>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2ALowModesSchurDiagOne(FMat &action);
    void makeLowModeV(FermionField &vout,
                      const FermionField &evec, const Real &eval);
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d,
                        const FermionField &evec, const Real &eval);
    void makeLowModeW(FermionField &wout,
                      const FermionField &evec, const Real &eval);
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d,
                        const FermionField &evec, const Real &eval);
    // See the note on A2ALowModesSchurDiagTwo::op().
    SchurOperatorBase<FermionField>& op(void);
private:
    using A2ALowModesSchurBase<FImpl>::action_;
    using A2ALowModesSchurBase<FImpl>::fGrid_;
    using A2ALowModesSchurBase<FImpl>::frbGrid_;
    using A2ALowModesSchurBase<FImpl>::tmp5_;
    using A2ALowModesSchurBase<FImpl>::src_o_;
    using A2ALowModesSchurBase<FImpl>::sol_e_;
    using A2ALowModesSchurBase<FImpl>::sol_o_;
    using A2ALowModesSchurBase<FImpl>::tmp_;
    SchurDiagOneOperator<FMat, FermionField> op_;
};

/******************************************************************************
 *                  Methods for V & W all-to-all vectors I/O                  *
 ******************************************************************************/
class A2AVectorsIo
{
public:
    struct Record: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(Record,
                                        unsigned int, index);
        Record(void): index(0) {}
    };
public:
    template <typename Field>
    static void write(const std::string fileStem, std::vector<Field> &vec,
                      const bool multiFile, const int trajectory = -1);
    // Write a single element under an explicit, caller-supplied index, without
    // ever holding a full std::vector<Field> resident -- used by modules that
    // stream their output bin by bin. Always writes one file per element
    // (i.e. the multiFile=true layout of write() above). The single-growing-
    // file (multiFile=false) counterpart is openWriter/writeRecord below: the
    // caller holds the ScidacWriter open across calls and streams records
    // into one file.
    template <typename Field>
    static void writeElement(const std::string fileStem, Field &elem,
                             const unsigned int index, const int trajectory = -1);
    template <typename Field>
    static void read(std::vector<Field> &vec, const std::string fileStem,
                     const bool multiFile, const int trajectory = -1);
    // Streaming access to the single-file (multiFile=false) layout, one
    // record at a time, so callers never hold a full std::vector<Field>:
    // open once, then write/read records in index order. readRecord is
    // sequential-only (SciDAC records cannot be seeked by index; the index
    // argument is a consistency check on the record read). readElement is the
    // one-element random-access read for the multiFile=true layout (the read
    // counterpart of writeElement above). All of these match the on-disk
    // layout of write()/read() exactly.
    static void openWriter(ScidacWriter &writer, const std::string fileStem,
                           GridBase *grid, const int trajectory = -1);
    template <typename Field>
    static void writeRecord(ScidacWriter &writer, Field &field,
                            const unsigned int index);
    static void openReader(ScidacReader &reader, const std::string fileStem,
                           const int trajectory = -1);
    template <typename Field>
    static void readRecord(ScidacReader &reader, Field &field,
                           const unsigned int index);
    template <typename Field>
    static void readElement(const std::string fileStem, Field &field,
                            const unsigned int index,
                            const int trajectory = -1);
private:
    static inline std::string vecFilename(const std::string stem, const int traj, 
                                          const bool multiFile)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        if (multiFile)
        {
            return stem + t;
        }
        else
        {
            return stem + t + ".bin";
        }
    }

    static inline std::string elementFilename(const std::string stem,
                                              const int traj,
                                              const unsigned int index)
    {
        return vecFilename(stem, traj, true) + "/elem"
               + std::to_string(index) + ".bin";
    }
};

/******************************************************************************
 *                    A2AHighModes template implementation                    *
 ******************************************************************************/
template <typename FImpl>
A2AHighModes<FImpl>::A2AHighModes(FMat &action, Solver &solver)
: action_(action)
, solver_(solver)
, fGrid_(action_.FermionGrid())
, tmp5_(fGrid_)
{}

template <typename FImpl>
void A2AHighModes<FImpl>::makeHighModeV(FermionField &vout,
                                        const FermionField &noise)
{
    solver_(vout, noise);
}

template <typename FImpl>
void A2AHighModes<FImpl>::makeHighModeV5D(FermionField &vout_4d,
                                          FermionField &vout_5d,
                                          const FermionField &noise)
{
    if (noise.Grid()->Dimensions() == fGrid_->Dimensions() - 1)
    {
        action_.ImportPhysicalFermionSource(noise, tmp5_);
    }
    else
    {
        tmp5_ = noise;
    }
    makeHighModeV(vout_5d, tmp5_);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2AHighModes<FImpl>::makeHighModeV(std::vector<FermionField> &vout,
                                        const std::vector<FermionField> &noise)
{
    if (vout.size() != noise.size())
    {
        HADRONS_ERROR(Size, "solution/noise batch size mismatch ("
                      + std::to_string(vout.size()) + " vs "
                      + std::to_string(noise.size()) + ")");
    }
    solver_(vout, noise);
}

template <typename FImpl>
void A2AHighModes<FImpl>::makeHighModeV5D(std::vector<FermionField> &vout_4d,
                                          std::vector<FermionField> &vout_5d,
                                          const std::vector<FermionField> &noise)
{
    unsigned int nBatch = noise.size();

    if ((vout_4d.size() != nBatch) || (vout_5d.size() != nBatch))
    {
        HADRONS_ERROR(Size, "solution/noise batch size mismatch");
    }

    std::vector<FermionField> src5(nBatch, fGrid_);

    for (unsigned int j = 0; j < nBatch; ++j)
    {
        if (noise[j].Grid()->Dimensions() == fGrid_->Dimensions() - 1)
        {
            action_.ImportPhysicalFermionSource(noise[j], src5[j]);
        }
        else
        {
            src5[j] = noise[j];
        }
    }
    makeHighModeV(vout_5d, src5);
    for (unsigned int j = 0; j < nBatch; ++j)
    {
        action_.ExportPhysicalFermionSolution(vout_5d[j], vout_4d[j]);
    }
}

template <typename FImpl>
void A2AHighModes<FImpl>::makeHighModeW(FermionField &wout,
                                        const FermionField &noise)
{
    wout = noise;
}

template <typename FImpl>
void A2AHighModes<FImpl>::makeHighModeW5D(FermionField &wout_4d,
                                          FermionField &wout_5d,
                                          const FermionField &noise)
{
    if (noise.Grid()->Dimensions() == fGrid_->Dimensions() - 1)
    {
        action_.ImportUnphysicalFermion(noise, wout_5d);
        wout_4d = noise;
    }
    else
    {
        wout_5d = noise;
        action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
    }
}

/******************************************************************************
 *                A2ALowModesSchurBase template implementation                *
 ******************************************************************************/
template <typename FImpl>
A2ALowModesSchurBase<FImpl>::A2ALowModesSchurBase(FMat &action)
: action_(action)
, fGrid_(action_.FermionGrid())
, frbGrid_(action_.FermionRedBlackGrid())
, tmp5_(fGrid_)
, src_o_(frbGrid_)
, sol_e_(frbGrid_)
, sol_o_(frbGrid_)
, tmp_(frbGrid_)
{}

/******************************************************************************
 *               A2ALowModesSchurDiagTwo template implementation               *
 ******************************************************************************/
template <typename FImpl>
A2ALowModesSchurDiagTwo<FImpl>::A2ALowModesSchurDiagTwo(FMat &action)
: A2ALowModesSchurBase<FImpl>(action)
, op_(action)
{}

template <typename FImpl>
void A2ALowModesSchurDiagTwo<FImpl>::makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval)
{
    src_o_ = evec;
    src_o_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_e_, vout);
    pickCheckerboard(Odd, sol_o_, vout);

    /////////////////////////////////////////////////////
    // v_ie = -(1/eval_i) * MeeInv Meo MooInv evec_i
    /////////////////////////////////////////////////////
    action_.MooeeInv(src_o_, tmp_);
    assert(tmp_.Checkerboard() == Odd);
    action_.Meooe(tmp_, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    action_.MooeeInv(sol_e_, tmp_);
    assert(tmp_.Checkerboard() == Even);
    sol_e_ = (-1.0 / eval) * tmp_;
    assert(sol_e_.Checkerboard() == Even);

    /////////////////////////////////////////////////////
    // v_io = (1/eval_i) * MooInv evec_i
    /////////////////////////////////////////////////////
    action_.MooeeInv(src_o_, tmp_);
    assert(tmp_.Checkerboard() == Odd);
    sol_o_ = (1.0 / eval) * tmp_;
    assert(sol_o_.Checkerboard() == Odd);
    setCheckerboard(vout, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    setCheckerboard(vout, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
}

template <typename FImpl>
void A2ALowModesSchurDiagTwo<FImpl>::makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &evec, const Real &eval)
{
    makeLowModeV(vout_5d, evec, eval);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2ALowModesSchurDiagTwo<FImpl>::makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval)
{
    src_o_ = evec;
    src_o_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_e_, wout);
    pickCheckerboard(Odd, sol_o_, wout);

    /////////////////////////////////////////////////////
    // w_ie = - MeeInvDag MoeDag Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_o_, tmp_);
    assert(tmp_.Checkerboard() == Odd);
    action_.MeooeDag(tmp_, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    action_.MooeeInvDag(sol_e_, tmp_);
    assert(tmp_.Checkerboard() == Even);
    sol_e_ = (-1.0) * tmp_;

    /////////////////////////////////////////////////////
    // w_io = Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_o_, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
    setCheckerboard(wout, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    setCheckerboard(wout, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
}

template <typename FImpl>
void A2ALowModesSchurDiagTwo<FImpl>::makeLowModeW5D(FermionField &wout_4d, 
                                                   FermionField &wout_5d, 
                                                   const FermionField &evec, 
                                                   const Real &eval)
{
    makeLowModeW(tmp5_, evec, eval);
    action_.DminusDag(tmp5_, wout_5d);
    action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
}

template <typename FImpl>
SchurOperatorBase<typename FImpl::FermionField>& A2ALowModesSchurDiagTwo<FImpl>::op(void)
{
    return op_;
}

/******************************************************************************
 *               A2ALowModesSchurDiagOne template implementation               *
 ******************************************************************************/
template <typename FImpl>
A2ALowModesSchurDiagOne<FImpl>::A2ALowModesSchurDiagOne(FMat &action)
: A2ALowModesSchurBase<FImpl>(action)
, op_(action)
{}

template <typename FImpl>
void A2ALowModesSchurDiagOne<FImpl>::makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval)
{
    src_o_ = evec;
    src_o_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_e_, vout);
    pickCheckerboard(Odd, sol_o_, vout);

    /////////////////////////////////////////////////////
    // v_ie = -(1/eval_i) * MeeInv Meo evec_i
    /////////////////////////////////////////////////////
    action_.Meooe(src_o_, tmp_);
    assert(tmp_.Checkerboard() == Even);
    action_.MooeeInv(tmp_, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    sol_e_ = (-1.0 / eval) * sol_e_;
    assert(sol_e_.Checkerboard() == Even);

    /////////////////////////////////////////////////////
    // v_io = (1/eval_i) * evec_i
    /////////////////////////////////////////////////////
    sol_o_ = (1.0 / eval) * src_o_;
    assert(sol_o_.Checkerboard() == Odd);
    setCheckerboard(vout, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    setCheckerboard(vout, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
}

template <typename FImpl>
void A2ALowModesSchurDiagOne<FImpl>::makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &evec, const Real &eval)
{
    makeLowModeV(vout_5d, evec, eval);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2ALowModesSchurDiagOne<FImpl>::makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval)
{
    src_o_ = evec;
    src_o_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_e_, wout);
    pickCheckerboard(Odd, sol_o_, wout);

    /////////////////////////////////////////////////////
    // w_io = MooInvDag Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_o_, tmp_);
    assert(tmp_.Checkerboard() == Odd);
    action_.MooeeInvDag(tmp_, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);

    /////////////////////////////////////////////////////
    // w_ie = - MeeInvDag MoeDag w_io
    /////////////////////////////////////////////////////
    action_.MeooeDag(sol_o_, tmp_);
    assert(tmp_.Checkerboard() == Even);
    action_.MooeeInvDag(tmp_, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    sol_e_ = (-1.0) * sol_e_;
    assert(sol_e_.Checkerboard() == Even);

    setCheckerboard(wout, sol_e_);
    assert(sol_e_.Checkerboard() == Even);
    setCheckerboard(wout, sol_o_);
    assert(sol_o_.Checkerboard() == Odd);
}

template <typename FImpl>
void A2ALowModesSchurDiagOne<FImpl>::makeLowModeW5D(FermionField &wout_4d,
                                                   FermionField &wout_5d,
                                                   const FermionField &evec,
                                                   const Real &eval)
{
    makeLowModeW(tmp5_, evec, eval);
    action_.DminusDag(tmp5_, wout_5d);
    action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
}

template <typename FImpl>
SchurOperatorBase<typename FImpl::FermionField>& A2ALowModesSchurDiagOne<FImpl>::op(void)
{
    return op_;
}

/******************************************************************************
 *               all-to-all vectors I/O template implementation               *
 ******************************************************************************/
template <typename Field>
void A2AVectorsIo::write(const std::string fileStem, std::vector<Field> &vec, 
                         const bool multiFile, const int trajectory)
{
    Record       record;
    GridBase     *grid = vec[0].Grid();
    ScidacWriter binWriter(grid->IsBoss());
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Writing vector " << i << std::endl;
            makeFileDir(fullFilename, grid);
            binWriter.open(fullFilename);
            record.index = i;
            binWriter.writeScidacFieldRecord(vec[i], record);
            binWriter.close();
        }
    }
    else
    {
        makeFileDir(filename, grid);
        binWriter.open(filename);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            LOG(Message) << "Writing vector " << i << std::endl;
            record.index = i;
            binWriter.writeScidacFieldRecord(vec[i], record);
        }
        binWriter.close();
    }
}

template <typename Field>
void A2AVectorsIo::writeElement(const std::string fileStem, Field &elem,
                                const unsigned int index, const int trajectory)
{
    GridBase     *grid = elem.Grid();
    ScidacWriter binWriter(grid->IsBoss());
    std::string  filename = elementFilename(fileStem, trajectory, index);

    makeFileDir(filename, grid);
    binWriter.open(filename);
    writeRecord(binWriter, elem, index);
    binWriter.close();
}

template <typename Field>
void A2AVectorsIo::read(std::vector<Field> &vec, const std::string fileStem,
                        const bool multiFile, const int trajectory)
{
    Record       record;
    ScidacReader binReader;
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Reading vector " << i << std::endl;
            binReader.open(fullFilename);
            binReader.readScidacFieldRecord(vec[i], record);
            binReader.close();
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
    }
    else
    {
        binReader.open(filename);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            LOG(Message) << "Reading vector " << i << std::endl;
            binReader.readScidacFieldRecord(vec[i], record);
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
        binReader.close();
    }
}

inline void A2AVectorsIo::openWriter(ScidacWriter &writer,
                                     const std::string fileStem,
                                     GridBase *grid,
                                     const int trajectory)
{
    std::string filename = vecFilename(fileStem, trajectory, false);

    makeFileDir(filename, grid);
    writer.open(filename);
}

template <typename Field>
void A2AVectorsIo::writeRecord(ScidacWriter &writer, Field &field,
                               const unsigned int index)
{
    Record record;

    LOG(Message) << "Writing vector " << index << std::endl;
    record.index = index;
    writer.writeScidacFieldRecord(field, record);
}

inline void A2AVectorsIo::openReader(ScidacReader &reader,
                                     const std::string fileStem,
                                     const int trajectory)
{
    std::string filename = vecFilename(fileStem, trajectory, false);

    reader.open(filename);
}

template <typename Field>
void A2AVectorsIo::readRecord(ScidacReader &reader, Field &field,
                              const unsigned int index)
{
    Record record;

    LOG(Message) << "Reading vector " << index << std::endl;
    reader.readScidacFieldRecord(field, record);
    if (record.index != index)
    {
        HADRONS_ERROR(Io, "vector index mismatch");
    }
}

template <typename Field>
void A2AVectorsIo::readElement(const std::string fileStem, Field &field,
                               const unsigned int index,
                               const int trajectory)
{
    ScidacReader reader;
    std::string  filename = elementFilename(fileStem, trajectory, index);

    reader.open(filename);
    readRecord(reader, field, index);
    reader.close();
}

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_
