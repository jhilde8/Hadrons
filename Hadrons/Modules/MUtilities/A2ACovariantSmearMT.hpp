#ifndef Hadrons_MUtilities_A2ACovariantSmearMT_hpp_
#define Hadrons_MUtilities_A2ACovariantSmearMT_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Grid/qcd/utils/CovariantSmearing.h>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       A2ACovariantSmearMT                                  *
 *  Reference implementation using GaussianSmear (sequential Cshift calls).   *
 *  Kept for regression testing against A2ACovariantSmear (stencil path).     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ACovariantSmearMTPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ACovariantSmearMTPar,
                                    std::string, a2aVectors,
                                    std::string, gauge,
                                    double, alpha,
                                    unsigned int, N,
                                    unsigned int, orthog_axis,
                                    std::string, output,
                                    bool, multiFile);
};

template <typename GImpl, typename FImpl>
class TA2ACovariantSmearMT: public Module<A2ACovariantSmearMTPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TA2ACovariantSmearMT(const std::string name);
    virtual ~TA2ACovariantSmearMT(void) {};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2ACovariantSmearMT, ARG(TA2ACovariantSmearMT<GIMPL, FIMPL>), MUtilities);

/******************************************************************************
 *                 TA2ACovariantSmearMT implementation                        *
 ******************************************************************************/
template <typename GImpl, typename FImpl>
TA2ACovariantSmearMT<GImpl, FImpl>::TA2ACovariantSmearMT(const std::string name)
: Module<A2ACovariantSmearMTPar>(name)
{}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmearMT<GImpl, FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().a2aVectors, par().gauge};

    return in;
}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmearMT<GImpl, FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename GImpl, typename FImpl>
void TA2ACovariantSmearMT<GImpl, FImpl>::setup(void)
{
    auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);
    envCreate(std::vector<FermionField>, getName(), 1, a2aVec.size(),
              envGetGrid(FermionField));
}

template <typename GImpl, typename FImpl>
void TA2ACovariantSmearMT<GImpl, FImpl>::execute(void)
{
    const auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);
    unsigned int Ni = a2aVec.size();

    const auto &U = envGet(GaugeField, par().gauge);

    std::vector<LatticeColourMatrix> Umu(Nd, env().getGrid());
    for (int i = 0; i < Nd; i++)
        Umu[i] = PeekIndex<LorentzIndex>(U, i);

    auto &a2aSmr = envGet(std::vector<FermionField>, getName());

    RealD width = par().alpha;
    unsigned int iterations = par().N;
    unsigned int orthog_axis = par().orthog_axis;

    LOG(Message) << "Starting A2A covariant smearing (MT/Cshift path), alpha=" << width
                 << " N=" << iterations << " orthog=" << orthog_axis << std::endl;
    LOG(Message) << par().a2aVectors << " has size " << Ni << std::endl;

    CovariantSmearing<GImpl> CovSmear;

    startTimer("GaussianSmear");
    for (int i = 0; i < (int)Ni; i++)
    {
        a2aSmr[i] = a2aVec[i];
        CovSmear.template GaussianSmear<FermionField>(Umu, a2aSmr[i], width, iterations, orthog_axis);
    }
    stopTimer("GaussianSmear");

    if (!par().output.empty())
    {
        startTimer("I/O");
        A2AVectorsIo::write(par().output, a2aSmr, par().multiFile, vm().getTrajectory());
        stopTimer("I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ACovariantSmearMT_hpp_
