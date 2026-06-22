#ifndef Hadrons_MUtilities_A2ACovariantSmear_hpp_
#define Hadrons_MUtilities_A2ACovariantSmear_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         A2ACovariantSmear                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

// Lifted from Grid/examples/Example_Laplacian.cc (CovariantLaplacianStencil).
// Single HaloExchange per application instead of one Cshift per direction.
#define A2ACOVSMEAR_LEG_LOAD(Dir)                                               \
    SE = st.GetEntry(ptype, Dir, ss);                                           \
    if (SE->_is_local) {                                                        \
        int perm = SE->_permute;                                                \
        chi = coalescedReadPermute(in[SE->_offset], ptype, perm, lane);         \
    } else {                                                                    \
        chi = coalescedRead(buf[SE->_offset], lane);                            \
    }                                                                           \
    acceleratorSynchronise();

#define A2ACOVSMEAR_LEG_LOAD_MULT(leg, polarisation)    \
    UU  = coalescedRead(U[ss](polarisation));            \
    A2ACOVSMEAR_LEG_LOAD(leg);                           \
    mult(&Uchi(), &UU, &chi());                          \
    res = res + Uchi;

template <typename GImpl, typename Field>
class CovariantLaplacianStencil
{
public:
    INHERIT_GIMPL_TYPES(GImpl);

    typedef typename Field::vector_object siteObject;

    template <typename vtype>
    using iImplDoubledGaugeField = iVector<iScalar<iMatrix<vtype, Nc>>, Nds>;
    typedef iImplDoubledGaugeField<vComplex> SiteDoubledGaugeField;
    typedef Lattice<SiteDoubledGaugeField> DoubledGaugeField;

    typedef CartesianStencil<siteObject, siteObject, SimpleStencilParams> StencilImpl;
    SimpleStencilParams p;

    GridBase *grid;
    StencilImpl Stencil;
    SimpleCompressor<siteObject> Compressor;
    DoubledGaugeField Uds;

    CovariantLaplacianStencil(GaugeField &Umu)
        : grid(Umu.Grid()),
          Stencil(grid, 6, Even,
                  std::vector<int>({Xdir, Ydir, Zdir, Xdir, Ydir, Zdir}),
                  std::vector<int>({1, 1, 1, -1, -1, -1}), p),
          Uds(grid)
    {
        for (int mu = 0; mu < Nd; mu++)
        {
            auto U = PeekIndex<LorentzIndex>(Umu, mu);
            PokeIndex<LorentzIndex>(Uds, U, mu);
            U = adj(Cshift(U, mu, -1));
            PokeIndex<LorentzIndex>(Uds, U, mu + 4);
        }
    }

    void M(const Field &_in, Field &_out)
    {
        Stencil.HaloExchange(_in, Compressor);

        auto st  = Stencil.View(AcceleratorRead);
        auto buf = st.CommBuf();
        autoView(in,  _in,  AcceleratorRead);
        autoView(out, _out, AcceleratorWrite);
        autoView(U,   Uds,  AcceleratorRead);

        typedef typename Field::vector_object        vobj;
        typedef decltype(coalescedRead(in[0]))       calcObj;
        typedef decltype(coalescedRead(U[0](0)))     calcLink;

        const int      Nsimd = vobj::Nsimd();
        const uint64_t NN    = grid->oSites();

        accelerator_for(ss, NN, Nsimd, {
            StencilEntry *SE;
            const int lane = acceleratorSIMTlane(Nsimd);
            calcObj  chi, res, Uchi;
            calcLink UU;
            int ptype;

            res = coalescedRead(in[ss]) * (-6.0);

            A2ACOVSMEAR_LEG_LOAD_MULT(0, Xp);
            A2ACOVSMEAR_LEG_LOAD_MULT(1, Yp);
            A2ACOVSMEAR_LEG_LOAD_MULT(2, Zp);
            A2ACOVSMEAR_LEG_LOAD_MULT(3, Xm);
            A2ACOVSMEAR_LEG_LOAD_MULT(4, Ym);
            A2ACOVSMEAR_LEG_LOAD_MULT(5, Zm);

            coalescedWrite(out[ss], res, lane);
        });
    }
};

#undef A2ACOVSMEAR_LEG_LOAD_MULT
#undef A2ACOVSMEAR_LEG_LOAD

class A2ACovariantSmearPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ACovariantSmearPar,
                                    std::string, a2aVectors,
                                    std::string, gauge,
                                    double, alpha,
                                    unsigned int, N,
                                    unsigned int, orthog_axis,
                                    std::string, output,
                                    bool, multiFile);
};

template <typename GImpl, typename FImpl>
class TA2ACovariantSmear: public Module<A2ACovariantSmearPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ACovariantSmear(const std::string name);
    // destructor
    virtual ~TA2ACovariantSmear(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2ACovariantSmear, ARG(TA2ACovariantSmear<GIMPL, FIMPL>), MUtilities);

/******************************************************************************
 *                 TA2ACovariantSmear implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
TA2ACovariantSmear<GImpl, FImpl>::TA2ACovariantSmear(const std::string name)
: Module<A2ACovariantSmearPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmear<GImpl, FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().a2aVectors, par().gauge};

    return in;
}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2ACovariantSmear<GImpl, FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2ACovariantSmear<GImpl, FImpl>::setup(void)
{
    auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);
    envCreate(std::vector<FermionField>, getName(), 1, a2aVec.size(),
              envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2ACovariantSmear<GImpl, FImpl>::execute(void)
{
    const auto &a2aVec = envGet(std::vector<FermionField>, par().a2aVectors);
    unsigned int Ni = a2aVec.size();

    auto &U = envGet(GaugeField, par().gauge);
    auto &a2aSmr = envGet(std::vector<FermionField>, getName());

    RealD width = par().alpha;
    unsigned int iterations = par().N;
    unsigned int orthog_axis = par().orthog_axis;
    Real coeff = (width * width) / Real(4 * iterations);

    LOG(Message) << "Starting A2A covariant smearing, alpha=" << width
                 << " N=" << iterations << " orthog=" << orthog_axis << std::endl;
    LOG(Message) << par().a2aVectors << " has size " << Ni << std::endl;

    CovariantLaplacianStencil<GImpl, FermionField> Lap(U);
    FermionField tmp(envGetGrid(FermionField));

    startTimer("GaussianSmear");
    for (int i = 0; i < (int)Ni; i++)
    {
        a2aSmr[i] = a2aVec[i];
        for (unsigned int n = 0; n < iterations; n++)
        {
            Lap.M(a2aSmr[i], tmp);
            a2aSmr[i] = a2aSmr[i] + coeff * tmp;
        }
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

#endif // Hadrons_MUtilities_A2ACovariantSmear_hpp_
