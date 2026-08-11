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
 *  Common interface for the Schur-convention-specific (DiagOne/DiagTwo)     *
 *  A2A vector constructions below. Lets a module hold a single              *
 *  std::unique_ptr<A2AVectorsSchurBase<FImpl>>, chosen and constructed at   *
 *  runtime from an XML schurConvention string in setup() (mirroring how    *
 *  MSolver::MixedPrecisionRBPrecCG and MSolver::A2AVectorsCoarseLow's       *
 *  checkOp_ already pick a concrete SchurDiagOneOperator/SchurDiagTwoOperator*
 *  at runtime), rather than needing the concrete type fixed at compile     *
 *  time. The high-mode methods never depend on the Schur convention (see   *
 *  A2AVectorsSchurDiagTwo/One below -- neither one's makeHighMode* touches *
 *  op_), so they are implemented once, here, rather than duplicated in     *
 *  both derived classes; only the low-mode methods are pure virtual.       *
 ******************************************************************************/
template <typename FImpl>
class A2AVectorsSchurBase
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2AVectorsSchurBase(FMat &action, Solver &solver);
    virtual ~A2AVectorsSchurBase(void) = default;
    virtual void makeLowModeV(FermionField &vout,
                              const FermionField &evec, const Real &eval) = 0;
    virtual void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d,
                                const FermionField &evec, const Real &eval) = 0;
    virtual void makeLowModeW(FermionField &wout,
                              const FermionField &evec, const Real &eval) = 0;
    virtual void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d,
                                const FermionField &evec, const Real &eval) = 0;
    // Exposes the concrete DiagOne/DiagTwo Schur operator used internally by
    // makeLowModeW/op_, so callers (e.g. A2AVectorsCoarseLow's in-program
    // eigenvector check) can validate evec_i against the exact same operator
    // this class uses to build V/W, rather than independently re-deriving
    // which SchurDiagOneOperator/SchurDiagTwoOperator to build from a second
    // copy of the schurConvention dispatch logic.
    virtual SchurOperatorBase<FermionField>& op(void) = 0;
    void makeHighModeV(FermionField &vout, const FermionField &noise);
    void makeHighModeV5D(FermionField &vout_4d, FermionField &vout_5d,
                         const FermionField &noise_5d);
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
 *                 Class to generate V & W all-to-all vectors                 *
 ******************************************************************************/
template <typename FImpl>
class A2AVectorsSchurDiagTwo : public A2AVectorsSchurBase<FImpl>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2AVectorsSchurDiagTwo(FMat &action, Solver &solver);
    ~A2AVectorsSchurDiagTwo(void) override = default;
    void makeLowModeV(FermionField &vout,
                      const FermionField &evec, const Real &eval) override;
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d,
                        const FermionField &evec, const Real &eval) override;
    void makeLowModeW(FermionField &wout,
                      const FermionField &evec, const Real &eval) override;
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d,
                        const FermionField &evec, const Real &eval) override;
    SchurOperatorBase<FermionField>& op(void) override;
private:
    using A2AVectorsSchurBase<FImpl>::action_;
    using A2AVectorsSchurBase<FImpl>::solver_;
    using A2AVectorsSchurBase<FImpl>::fGrid_;
    using A2AVectorsSchurBase<FImpl>::tmp5_;
    GridBase                                 *frbGrid_, *gGrid_;
    bool                                     is5d_;
    FermionField                             src_o_, sol_e_, sol_o_, tmp_;
    SchurDiagTwoOperator<FMat, FermionField> op_;
};

template <typename FImpl>
class A2AVectorsSchurDiagOne : public A2AVectorsSchurBase<FImpl>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    A2AVectorsSchurDiagOne(FMat &action, Solver &solver);
    ~A2AVectorsSchurDiagOne(void) override = default;
    void makeLowModeV(FermionField &vout,
                      const FermionField &evec, const Real &eval) override;
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d,
                        const FermionField &evec, const Real &eval) override;
    void makeLowModeW(FermionField &wout,
                      const FermionField &evec, const Real &eval) override;
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d,
                        const FermionField &evec, const Real &eval) override;
    SchurOperatorBase<FermionField>& op(void) override;
private:
    using A2AVectorsSchurBase<FImpl>::action_;
    using A2AVectorsSchurBase<FImpl>::solver_;
    using A2AVectorsSchurBase<FImpl>::fGrid_;
    using A2AVectorsSchurBase<FImpl>::tmp5_;
    GridBase                                 *frbGrid_, *gGrid_;
    bool                                     is5d_;
    FermionField                             src_o_, sol_e_, sol_o_, tmp_;
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
    // (i.e. the multiFile=true layout of write() above); there is no
    // single-growing-file (multiFile=false) counterpart, since that would
    // require keeping a ScidacWriter open across calls.
    template <typename Field>
    static void writeElement(const std::string fileStem, Field &elem,
                             const unsigned int index, const int trajectory = -1);
    template <typename Field>
    static void read(std::vector<Field> &vec, const std::string fileStem,
                     const bool multiFile, const int trajectory = -1);
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
};

/******************************************************************************
 *                 A2AVectorsSchurBase template implementation                *
 ******************************************************************************/
template <typename FImpl>
A2AVectorsSchurBase<FImpl>::A2AVectorsSchurBase(FMat &action, Solver &solver)
: action_(action)
, solver_(solver)
, fGrid_(action_.FermionGrid())
, tmp5_(fGrid_)
{}

template <typename FImpl>
void A2AVectorsSchurBase<FImpl>::makeHighModeV(FermionField &vout,
                                               const FermionField &noise)
{
    solver_(vout, noise);
}

template <typename FImpl>
void A2AVectorsSchurBase<FImpl>::makeHighModeV5D(FermionField &vout_4d,
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
void A2AVectorsSchurBase<FImpl>::makeHighModeW(FermionField &wout,
                                               const FermionField &noise)
{
    wout = noise;
}

template <typename FImpl>
void A2AVectorsSchurBase<FImpl>::makeHighModeW5D(FermionField &wout_4d,
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
 *               A2AVectorsSchurDiagTwo template implementation               *
 ******************************************************************************/
template <typename FImpl>
A2AVectorsSchurDiagTwo<FImpl>::A2AVectorsSchurDiagTwo(FMat &action, Solver &solver)
: A2AVectorsSchurBase<FImpl>(action, solver)
, frbGrid_(action.FermionRedBlackGrid())
, gGrid_(action.GaugeGrid())
, src_o_(frbGrid_)
, sol_e_(frbGrid_)
, sol_o_(frbGrid_)
, tmp_(frbGrid_)
, op_(action)
{}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval)
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
void A2AVectorsSchurDiagTwo<FImpl>::makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &evec, const Real &eval)
{
    makeLowModeV(vout_5d, evec, eval);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2AVectorsSchurDiagTwo<FImpl>::makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval)
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
void A2AVectorsSchurDiagTwo<FImpl>::makeLowModeW5D(FermionField &wout_4d, 
                                                   FermionField &wout_5d, 
                                                   const FermionField &evec, 
                                                   const Real &eval)
{
    makeLowModeW(tmp5_, evec, eval);
    action_.DminusDag(tmp5_, wout_5d);
    action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
}

template <typename FImpl>
SchurOperatorBase<typename FImpl::FermionField>& A2AVectorsSchurDiagTwo<FImpl>::op(void)
{
    return op_;
}

/******************************************************************************
 *               A2AVectorsSchurDiagOne template implementation               *
 ******************************************************************************/
template <typename FImpl>
A2AVectorsSchurDiagOne<FImpl>::A2AVectorsSchurDiagOne(FMat &action, Solver &solver)
: A2AVectorsSchurBase<FImpl>(action, solver)
, frbGrid_(action.FermionRedBlackGrid())
, gGrid_(action.GaugeGrid())
, src_o_(frbGrid_)
, sol_e_(frbGrid_)
, sol_o_(frbGrid_)
, tmp_(frbGrid_)
, op_(action)
{}

template <typename FImpl>
void A2AVectorsSchurDiagOne<FImpl>::makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval)
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
void A2AVectorsSchurDiagOne<FImpl>::makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &evec, const Real &eval)
{
    makeLowModeV(vout_5d, evec, eval);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2AVectorsSchurDiagOne<FImpl>::makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval)
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
void A2AVectorsSchurDiagOne<FImpl>::makeLowModeW5D(FermionField &wout_4d,
                                                   FermionField &wout_5d,
                                                   const FermionField &evec,
                                                   const Real &eval)
{
    makeLowModeW(tmp5_, evec, eval);
    action_.DminusDag(tmp5_, wout_5d);
    action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
}

template <typename FImpl>
SchurOperatorBase<typename FImpl::FermionField>& A2AVectorsSchurDiagOne<FImpl>::op(void)
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
    Record       record;
    GridBase     *grid = elem.Grid();
    ScidacWriter binWriter(grid->IsBoss());
    std::string  filename     = vecFilename(fileStem, trajectory, true);
    std::string  fullFilename = filename + "/elem" + std::to_string(index) + ".bin";

    LOG(Message) << "Writing vector " << index << std::endl;
    makeFileDir(fullFilename, grid);
    binWriter.open(fullFilename);
    record.index = index;
    binWriter.writeScidacFieldRecord(elem, record);
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

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_
