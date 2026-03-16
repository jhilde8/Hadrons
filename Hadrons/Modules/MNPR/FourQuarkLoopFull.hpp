/*
 * FourQuarkLoopFull.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Fabian Joswig <fabian.joswig@wwu.de>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Ryan Abbott <rabbott@mit.edu>
 * Author: Simon Bürger <simon.buerger@rwth-aachen.de>
 * Author: rabbott <rabbott4927@gmail.com>
 * Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
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
#ifndef Hadrons_MNPR_FourQuarkLoopFull_hpp_
#define Hadrons_MNPR_FourQuarkLoopFull_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MNPR/NPRUtils.hpp>


/******************************************************************************
 *                               FourQuarkLoop                                *
 ******************************************************************************/
BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MNPR)

/* This module computes loop diagrams for 4-quark operators for use in
 * non-perturbative renormalization. It can compute two diagrams, referred to
 * as "connected loop" and "disconnected loop". Diagramatically connected loop
 * refers to the following diagram:
 *
 *    qOut,pOut            qOut,pOut
 *        \                   /
 *         ^     _          /
 *          \  /   \      /
 *  \Gamma_A x      |    |
 *                  ^    ^
 *  \Gamma_B x      |    |
 *          /  \ _ /      \
 *         ^                \
 *        /                   \
 *    qIn,pIn              qIn,pIn
 *
 * Meanwhile the disconnected loop diagram is:
 *
 *     qOut,pOut                       qOut,pOut
 *         \                              /
 *          ^         _ _ _ _           /
 *           \      /         \       /
 *            \    |           |     |
 *   \Gamma_A  x   x \Gamma_B  ^     ^
 *            /    |           |     |
 *           /      \ _ _ _ _ /       \
 *          ^                           \
 *         /                              \
 *     qIn,pIn                          qIn,pIn
 *
 * In the diagrams above both 'x's are located at the same position for a given
 * diagram. At each 'x' the \Gamma_A and \Gamma_B are insertions of
 * gamma-matrix structures. There are many possiblities for the choice(s) of
 * \Gamma_A and \Gamma_B; the parameter gamma_basis controls which are used
 * (see comments below on gamma_basis).
 *
 * This diagram does *not* compute the spectator quarks (the lines on the right
 * above). The spectator quarks may be computed independently using the
 * ExternalLeg module.
 *
 * All diagrams are phased so that they carry 0 net momentum; consequently this
 * module is suitable for renormalizing 4-quark operators in the RI/SMOM scheme
 * where pIn != pOut.
 *
 * Outputs:
 *   There are two primary outputs of the module, corresponding to the dieagams with 2 external legs
 *   and 4 external legs, respoectively. The 2 quark external state is written to a file as vectors
 *   of SpinColourMatrix, the 4 quark external state is written to a file as a vector of
 *   SpinColourSpinColour matrix. In the four quark case, we can utilize a long distance truncation
 *   that uses a window function to suppress the noisiest contributions of these disconnected diagrams.
 *
 * Inputs:
 *   loop_type - should be either "connected" or "disconnected". Determines
 *               whether the module computes the connected or disconnected loop
 *               diagram listed above
 *
 *   qIn, qOut - Input and output propagators
 *
 *   pIn, pOut - The momenta corresponding to qIn and qOut, respectively
 *
 *   qLoop - A propagagtor field for the loop, presumably calculated using
 *              all-to-all vectors. One possible way to generate this input is
 *              using the A2ALoop module.
 *
 *   num_a2a_vectors - The number of all-to-all vectors used to generate
 *                     'qLoop'. This is necessary when using the module
 *                     A2ALoop since A2ALoop doesn't know the number of vectors
 *                     used, and consequently the loop calculated there is not
 *                     an average but rather a sum. We correct for this by
 *                     dividing the loop by num_a2a_vectors.
 *   TODO: See comment in the code, not sure whether we need/want this
 *
 *   truncation - boolean value that sets whether or not to do this truncation on the spectator. 
 *
 *   gamma_basis - determines the insertions \Gamma_A and \Gamma_B; see above
 *                 for more detail and see below for what options are available
 *
 *   output - string which determines where the outputs are placed
 *
 */


class FourQuarkLoopFullPar : Serializable
{
public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FourQuarkLoopFullPar,
                                        std::string, loop_type,
                                        std::string, qIn,
                                        std::string, qOut,
                                        std::string, pIn,
                                        std::string, pOut,
                                        std::string, qLoop,
                                        int, num_a2a_vectors,
					std::string, window,
					double, fwhm,
                                        std::string, gamma_basis,
                                        std::string, output);
};


template <typename FImpl>
class TFourQuarkLoopFull : public Module<FourQuarkLoopFullPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,)
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::vector<Gamma::Algebra>, gammaA,
                                        std::vector<Gamma::Algebra>, gammaB,
                                        std::string,  pIn,
                                        std::string,  pOut);
    };
    //TODO: I use correlator here, but the result is a SpinColourMatrix per gamma matrix, not per timeslice. Do we want another struct?
    typedef Correlator<Metadata, SpinColourMatrix> Result2q;
    typedef Correlator<Metadata, SpinColourSpinColourMatrix> Result4q; //full un-truncated result

    TFourQuarkLoopFull(const std::string name);
    virtual ~TFourQuarkLoopFull(void) {};

    virtual std::vector<std::string> getInput();
    virtual std::vector<std::string> getOutput();
    virtual std::vector<std::string> getOutputFiles(void);

protected:
    virtual void setup(void);
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FourQuarkLoopFull, ARG(TFourQuarkLoopFull<FIMPL>), MNPR);

template <typename FImpl>
TFourQuarkLoopFull<FImpl>::TFourQuarkLoopFull(const std::string name)
    : Module<FourQuarkLoopFullPar>(name)
{}

template <typename FImpl>
std::vector<std::string> TFourQuarkLoopFull<FImpl>::getInput()
{
    std::vector<std::string> in = { par().qIn, par().qOut, par().qLoop };

    return in;
}

template <typename FImpl>
std::vector<std::string> TFourQuarkLoopFull<FImpl>::getOutput()
{
    std::vector<std::string> out = {getName(), getName() + "_fourQuark", getName() + "_twoQuark", getName() + "_fourQuarkTrunc"};

    return out;
}

template <typename FImpl>
std::vector<std::string> TFourQuarkLoopFull<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};

    return output;
}

template <typename FImpl>
void TFourQuarkLoopFull<FImpl>::setup()
{
    LOG(Message) << "Running setup for FourQuarkLoopFull"
        << std::endl;
    if (par().loop_type != "connected" && par().loop_type != "disconnected")
    {
        HADRONS_ERROR(Definition, "Unkown loop type");
    }

    envTmpLat(ComplexField, "pDotXin");
    envTmpLat(ComplexField, "w_p");

    envTmpLat(PropagatorField, "loop");
    envTmpLat(PropagatorField, "bilinear");
    envTmpLat(PropagatorField, "spectator");
    envTmpLat(PropagatorField, "STilde"); 
    envTmpLat(PropagatorField, "wind_spec"); 
    envTmpLat(PropagatorField, "bilinear_sum");

    envTmpLat(ComplexField, "bilinear_phase");
    if (par().loop_type == "disconnected")
    {
        envTmpLat(ComplexField, "loop_trace");
    }

    envCreate(HadronsSerializable, getName() + "_fourQuark", 1, 0);
    envCreate(HadronsSerializable, getName() + "_twoQuark", 1, 0);
}

template <typename FImpl>
void TFourQuarkLoopFull<FImpl>::execute()
{
    LOG(Message) << "Computing contractions '" << getName()
        << "' using source propagators '" << par().qIn << "' and '" << par().qOut << "'"
        << std::endl;

    PropagatorField &qIn = envGet(PropagatorField, par().qIn);
    PropagatorField &qOut = envGet(PropagatorField, par().qOut);

    PropagatorField &loop_raw = envGet(PropagatorField, par().qLoop);
    envGetTmp(PropagatorField, loop);
    // The loop is supposed to be an average, but the A2ALoop module doesn't
    // divide out by the number of all-to-all vectors; therefore we must divide
    // here.
    // TODO: We don't really want this, do we? We probably don't even want to assume
    // A2ALoop as the only module which can produce those, and the normalisation can be
    // done later, or directly in the A2ALoop module?
    loop = loop_raw * (1.0 / par().num_a2a_vectors);

    envGetTmp(PropagatorField, bilinear); //usual bilineari
    envGetTmp(PropagatorField, bilinear_sum); //tmp prop field for tensor product (can be reused)
    envGetTmp(PropagatorField, spectator); //usual spectator 
    envGetTmp(PropagatorField, STilde); //holds FFT of spectator.
    //envGetTmp(PropagatorField, wind_spec); //holds iFFT of window * spectator. 

    Coordinate                  latt_size = GridDefaultLatt();
    std::vector<Real> pIn = strToVec<Real>(par().pIn);
    std::vector<Real> pOut = strToVec<Real>(par().pOut);

    Result4q fourq_result;
    Result2q twoq_result;
    fourq_result.info.pIn  = par().pIn;
    fourq_result.info.pOut = par().pOut;
    twoq_result.info.pIn  = par().pIn;
    twoq_result.info.pOut = par().pOut;

    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);

    envGetTmp(ComplexField, bilinear_phase);
    envGetTmp(ComplexField, pDotXin);
    envGetTmp(ComplexField, w_p);

    GridBase *grid = pDotXin.Grid(); //get the grid directly from env please. 

    Complex Ci(0.0, 1.0);

    bool trunc = false; //necessary for tensor product demarcation in compute_diagrams
    double fwhm = par().fwhm; 

    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++)
    {
        volume *= latt_size[mu];
    }

    //lambda function for the case of using a window function. Computes the FFT, convolution in momentum space, and iFFT. 
    auto spectator_window_convolution = [&]() {
    
        LOG(Message) << "Computing window-spectator convolution" << std::endl;

        FFT fft((GridCartesian *)grid);
        fft.FFT_all_dim(STilde, spectator, FFT::forward);
        STilde = w_p * STilde;
        fft.FFT_all_dim(spectator, STilde, FFT::backward); 
    
        //normalize the result by V/sum(w_p)
        Complex window_norm = sum(w_p); 
        spectator *= volume/window_norm; 
    };

    LOG(Message) << "Calculating phases" << std::endl;

    NPRUtils<FImpl>::phase(bilinear_phase,pIn,pOut);
    NPRUtils<FImpl>::dot(pDotXin, pOut); 
    spectator = qIn * exp(-Ci * pDotXin); //position space spectator field
    SpinColourMatrix s1 = sum(spectator); //this completes the momentum projection for the un-windowed spectator. 

    LOG(Message) << "Computing window functions" << std::endl;

    WindowType window = parseWindowType(par().window);
    switch(window) {
        case WindowType::RECTANGLE: 
	    trunc = true;
	    w_p = NPRUtils<FImpl>::getRectWindow(grid);
	    spectator_window_convolution();
	    break;
	case WindowType::GAUSSIAN:
	    trunc = true;
	    w_p = NPRUtils<FImpl>::getGaussianWindow(grid, fwhm);
    	    spectator_window_convolution();
	    break;
	case WindowType::NONE: //spectator is already ready, so we do nothing here. 
	default:
	    break; 
    }

    /*
    LOG(Message) << "Computing window-spectator convolution" << std::endl;

    //Here we actually use the window, so depending on the 'window' parameter we set one of the cached windows 
    //to the used w_p variable
	
    FFT fft((GridCartesian *)grid);
    fft.FFT_all_dim(STilde, spectator, FFT::forward);
    STilde = w_p * STilde;
    
    fft.FFT_all_dim(wind_spec, STilde, FFT::backward); 
    
    //normalize the result by the sum of the spectator
    Complex window_norm = sum(w_p); 
    wind_spec *= volume/window_norm; 
    */

    LOG(Message) << "Computing diagrams" << std::endl;

    const bool loop_disconnected = par().loop_type == "disconnected";

    auto compute_diagrams = [&](Gamma gamma_A, Gamma gamma_B, bool print = true) {

        twoq_result.info.gammaA.push_back(gamma_A.g);
        twoq_result.info.gammaB.push_back(gamma_B.g);
        fourq_result.info.gammaA.push_back(gamma_A.g);
        fourq_result.info.gammaB.push_back(gamma_B.g);

        if (print)
	{
            LOG(Message) << "Computing diagrams with GammaA = "
                << gamma_A.g << ", " << "GammaB = " << gamma_B.g
                << std::endl;
        }

        if (loop_disconnected)
	{
            // Disconnected loop diagram
            envGetTmp(ComplexField, loop_trace);
            loop_trace = trace(loop * gamma_B);
            bilinear = g5 * adj(qOut) * g5 * gamma_A * qIn;
            bilinear = (bilinear_phase * loop_trace) * bilinear;
        }
        else
	{
            // Connected loop diagram
            bilinear = bilinear_phase * (g5 * adj(qOut) * g5 * gamma_A * loop * gamma_B * qIn);
        }
        twoq_result.corr.push_back( (1.0 / volume) * sum(bilinear) );

         /*The quantity we actually need out of this is
         *
         * \sum_x bilinear_x \otimes spectator exp(-2i (pIn - pOut) \cdot x)
         *
         * Which is the appropriately phased version of the loop diagram so
         * that diagram carries zero net momentum. At this point we must apply one more 
	 * phase factor, then complete the position space sum, and the do the tensor product with the summed spectator
         *
         */
        bilinear = bilinear * bilinear_phase; //position space bilinear
	
	/*The truncation instroduces a window function to the spectator, which allows us to truncate the spectator contribution the further it is from the loop
	 * \sum_x b(x) \otimes \sum_y w(|x-y|) s(y)
	 * here we must take the tensor product here, since the spectator position sum depends on the location of the weak operator insertion. 
	 * This can be done efficiently with the FFT since this is a convolution, we do the FFT on the window and the spectator separately, 
	 * then do a product between the two in momentum space, iFFT the result, and complete the tensor product at each site,
	 *
	 * \sum_x b(x) \otimes g(x) , where g(x) = iFFT(w(p) * s(p)). 
	 * choice of window function is set in the input xml. This sets the trunc boolean. 
	 */
	if (trunc == true){
	    //truncated spectator-bilinear tensor product in position space. This uses the iFFT of the windowed spectator, and computes \sum_x b(x) \otimes iFFT(w * g)(x)
	    LOG(Message) << "Computing tensor products with windowed spectator" << std::endl;
	    SpinColourSpinColourMatrix lret;
	    bilinear_sum = Zero();
	    lret = NPRUtils<FImpl>::tensorProdSum(bilinear_sum, bilinear, spectator); 
	    fourq_result.corr.push_back((1.0 / volume) * lret);
	} else {
	    LOG(Message) << "Computing single tensor product with un-windowed spectator" << std::endl;
	    SpinColourMatrix b1 = sum(bilinear); //momentum space bilnear term

	    SpinColourSpinColourMatrix mom_ret;
	    NPRUtils<FImpl>::tensorSiteProd(mom_ret, b1, s1); //single site tensor product between momentum space spectator and bilinear

            fourq_result.corr.push_back( (1.0 / volume) * mom_ret );
	}

    };

    /* Choices of gamma_basis
     *
     * There are currently 4 choices for gamma_basis
     *
     * all - Computes every possible structure, 256 in total
     *
     * diagonal - Computes every structure where \Gamma_A = \Gamma_B, leading
     *            to 16 structures
     *
     * diagonal_va - Computes the diagrams needed to create vector/axial
     *               structure, i.e. structures of the form
     *               (qbar \Gamma^\mu q') (qbar'' \Gamma_\mu q''') where
     *               \Gamma^\mu = \gamma^\mu or \gamma^\mu \gamma^5.
     *
     * diagona_va_sp - Computes everything in 'diagonal_va' plus diagrams with
     *                 scalar/pseudoscalar strctures, i.e. structures of the
     *                 form (qbar \Gamma q') (qbar'' \Gamma q''') where \Gamma
     *                 = 1 or \gamma^5
     *
     * diagona_va_sp_tt - Computes everything in 'diagonal_va_sp' plus diagrams with
     *                    tensor strctures
     */

    std::string gamma_basis = par().gamma_basis;
    if (gamma_basis == "all")
    {
        for (Gamma gammaA: Gamma::gall)
	{
            for (Gamma gammaB: Gamma::gall)
	    {
                compute_diagrams(gammaA, gammaB);
            }
        }
    }
    else if (gamma_basis == "diagonal")
    {
        for (Gamma g: Gamma::gall)
	{
            compute_diagrams(g, g);
        }
    }
    else if (gamma_basis == "diagonal_va" || gamma_basis == "diagonal_va_sp" || gamma_basis == "diagonal_va_sp_tt")
    {
        for (int mu = 0; mu < 4; mu++)
	{
            Gamma gmu = Gamma::gmu[mu];
            Gamma gmug5 = Gamma::mul[gmu.g][Gamma::Algebra::Gamma5];
            compute_diagrams(gmu, gmu);
            compute_diagrams(gmu, gmug5);
            compute_diagrams(gmug5, gmu);
            compute_diagrams(gmug5, gmug5);
        }
        if (gamma_basis == "diagonal_va_sp" || gamma_basis == "diagonal_va_sp_tt")
	{
            Gamma identity = Gamma(Gamma::Algebra::Identity);

            compute_diagrams(identity, identity);
            compute_diagrams(identity, g5);
            compute_diagrams(g5, identity);
            compute_diagrams(g5, g5);
        }
	if (gamma_basis == "diagonal_va_sp_tt") {
            const std::array<const Gamma, 6> gsigma = {{
                  Gamma(Gamma::Algebra::SigmaXT),
                  Gamma(Gamma::Algebra::SigmaXY),
                  Gamma(Gamma::Algebra::SigmaXZ),
                  Gamma(Gamma::Algebra::SigmaYT),
                  Gamma(Gamma::Algebra::SigmaYZ),
                  Gamma(Gamma::Algebra::SigmaZT)}};

            for (Gamma gammaA: gsigma) {
                    compute_diagrams(gammaA, gammaA);
            }

            compute_diagrams(Gamma(Gamma::Algebra::SigmaXT), Gamma(Gamma::Algebra::SigmaYZ));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaXY), Gamma(Gamma::Algebra::SigmaZT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaXZ), Gamma(Gamma::Algebra::SigmaYT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaYT), Gamma(Gamma::Algebra::SigmaXZ));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaYZ), Gamma(Gamma::Algebra::SigmaXT));
            compute_diagrams(Gamma(Gamma::Algebra::SigmaZT), Gamma(Gamma::Algebra::SigmaXY));
        }
    }
    else
    {
        LOG(Error) << "Error: unkown gamma_basis: '" << gamma_basis << "'"
            << std::endl;
    }

    LOG(Message) << "Done computing loop diagrams" << std::endl;
    saveResult(par().output + "_fourQuark", "FourQuarkLoop", fourq_result);
    saveResult(par().output + "_twoQuark", "TwoQuarkLoop", twoq_result); 
    envGet(HadronsSerializable, getName() + "_fourQuark") = fourq_result;
    envGet(HadronsSerializable, getName() + "_twoQuark") = twoq_result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
