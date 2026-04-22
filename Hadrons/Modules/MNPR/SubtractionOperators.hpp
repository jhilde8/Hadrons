/*
 * SubtractionOperators.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
 * Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu> (windowing)
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
#ifndef Hadrons_MNPR_SubtractionOperators_hpp_
#define Hadrons_MNPR_SubtractionOperators_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MNPR/NPRUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SubtractionOperators                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNPR)

class SubtractionOperatorsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SubtractionOperatorsPar,
                                    std::string, qIn,
                                    std::string, qOut,
                                    std::string, pIn,
                                    std::string, pOut,
                                    std::string, gauge,
				    std::string, window,
				    double, fwhm,
                                    std::string, output);
};

template <typename FImpl>
class TSubtractionOperators: public Module<SubtractionOperatorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
        public:
            // Contains information on both the 2 and 4-quark external state
            // diagrams with the given subtraction operator
            class OperatorResult : Serializable {
                public:
                GRID_SERIALIZABLE_CLASS_MEMBERS(OperatorResult,
                        SpinColourSpinColourMatrix, fourq,
                        SpinColourMatrix, twoq);
            };
            GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                    OperatorResult, dslash_left,
                    OperatorResult, dslash_gamma5_left,
                    OperatorResult, dslash_right,
                    OperatorResult, dslash_gamma5_right,
                    OperatorResult, scalar,
                    OperatorResult, psuedoscalar);
    };
    // constructor
    TSubtractionOperators(const std::string name);
    // destructor
    virtual ~TSubtractionOperators(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);

    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SubtractionOperators, TSubtractionOperators<FIMPL>, MNPR);

/******************************************************************************
 *                 TSubtractionOperators implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSubtractionOperators<FImpl>::TSubtractionOperators(const std::string name)
: Module<SubtractionOperatorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSubtractionOperators<FImpl>::getInput(void)
{
    std::vector<std::string> in = { par().qOut, par().qIn, par().gauge };

    return in;
}

template <typename FImpl>
std::vector<std::string> TSubtractionOperators<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TSubtractionOperators<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};

    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSubtractionOperators<FImpl>::setup(void)
{
    LOG(Message) << "Running setup for SubtractionOperators" << std::endl;

    envTmpLat(PropagatorField, "Dslash_qIn");
    envTmpLat(PropagatorField, "Dslash_qOut");

    envTmpLat(PropagatorField, "bilinear");
    envTmpLat(PropagatorField, "bilinear_tmp");
    envTmpLat(PropagatorField, "spectator");
    envTmpLat(PropagatorField, "stilde"); //FFT of spectator

    envTmpLat(ComplexField, "bilinear_phase");
    envTmpLat(ComplexField, "pDotXOut");
    envTmpLat(ComplexField, "coordinate");
    envTmpLat(ComplexField, "w_p");

    envTmpLat(PropagatorField, "tmp");

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSubtractionOperators<FImpl>::execute(void)
{
    LOG(Message) << "Computing subtraction operators '"
        << "' using source propagators '" << par().qIn << "' and '" << par().qOut << "'"
        << std::endl;

    auto &Umu = envGet(GaugeField, par().gauge);

    PropagatorField qIn = envGet(PropagatorField, par().qIn);
    PropagatorField qOut = envGet(PropagatorField, par().qOut);

    envGetTmp(PropagatorField, Dslash_qIn);
    envGetTmp(PropagatorField, Dslash_qOut);

    envGetTmp(PropagatorField, bilinear);
    envGetTmp(PropagatorField, bilinear_tmp);
    envGetTmp(PropagatorField, spectator);
    envGetTmp(PropagatorField, stilde);

    Coordinate                  latt_size = GridDefaultLatt();
    std::vector<Real> pIn = strToVec<Real>(par().pIn);
    std::vector<Real> pOut = strToVec<Real>(par().pOut);

    envGetTmp(ComplexField, bilinear_phase);
    envGetTmp(ComplexField, pDotXOut);
    envGetTmp(ComplexField, coordinate);
    envGetTmp(ComplexField, w_p);

    Result result;
    Gamma g5(Gamma::Algebra::Gamma5);

    //GridBase *grid = getGrid5d(const bool redBlack = false, const unsigned int Ls = 16); 
    GridBase *grid = pDotXOut.Grid();

    bool trunc = false; //windowing flag
    double fwhm = par().fwhm; //this parameter is passed in regardless of if we window or not, we ignore it if trunc=false. 

    //// Compute volume
    Real volume = 1.0;
    for (int mu = 0; mu < Nd; mu++)
    {
        volume *= latt_size[mu];
    }

    LOG(Message) << "Calculating phases" << std::endl;

    NPRUtils<FImpl>::phase(bilinear_phase,pIn,pOut);
    NPRUtils<FImpl>::dot(pDotXOut,pOut);

    //// Compute Dslash for both propagators
    NPRUtils<FImpl>::dslash(Dslash_qIn, qIn, Umu);
    NPRUtils<FImpl>::dslash(Dslash_qOut, qOut, Umu);

    //// Compute spectator quark for 4-quark diagrams
    Complex Ci = Complex(0.0, 1.0);
    spectator = qIn * exp(-Ci * pDotXOut);
    SpinColourMatrixScalar spectator_sum = sum(spectator);

    //// Computing window function and spectator window convolution 
    LOG(Message) << "Computing window functions" << std::endl;

    auto spectator_window_convolution = [&]() {

        LOG(Message) << "Computing window-spectator convolution" << std::endl;

        FFT fft((GridCartesian *)grid);
        fft.FFT_all_dim(stilde, spectator, FFT::forward);
        stilde = w_p * stilde;
        fft.FFT_all_dim(spectator, stilde, FFT::backward); //spectator is now windowed 

        //normalize the result by V/sum(w_p)
        Complex window_norm = sum(w_p);
        spectator *= volume/window_norm;
    };

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
	case WindowType::DELTA:
	    trunc=true; 
	    //in the delta window case we don't need to waste time with a convolution since the position space 
	    //window has 1 at one site and zero everywhere else. Thus we just do the tensor product sum with the 
	    //position space spectator since the delta function contracts the sums.  
	    break;
	case WindowType::ZERO:
	    trunc=true;
	    spectator=Zero();
	    break;
	case WindowType::NONE: //spectator is already ready, so we do nothing here. 
	default:
	    break; 
    }

    //// Compute results
    auto compute_result = [&] (typename Result::OperatorResult &res)
    {

        bilinear = bilinear_phase * bilinear;
        res.twoq = (1.0 / volume) * sum(bilinear); //twoq result is same with and without window
        bilinear = bilinear_phase * bilinear; //second phase applied for fourq result
        
	if (trunc==true){
	    LOG(Message) << "Computing subtraction diagram with windowed spectator" << std::endl;	
	    //we want to use 'spectator' and 'bilinear' to do the tensor product at each site. 
	    bilinear_tmp = Zero(); //reset tmp prop field for tensor product functionality
	    res.fourq = (1.0 / volume) * NPRUtils<FImpl>::tensorProdSum(bilinear_tmp, bilinear, spectator);

	} else {
	    LOG(Message) << "Computing standard subtraction diagram" << std::endl;

            SpinColourMatrixScalar bilinear_avg = (1.0 / volume) * sum(bilinear);
            NPRUtils<FImpl>::tensorSiteProd(res.fourq, bilinear_avg, spectator_sum); //tensor product of summed propagator fields
	}

    };

    // The expression we want to compute here is
    //
    //  gamma^5 * adj(D_mu qOut) gamma^5 gamma^mu qIn
    //
    // But we first anti-commute the gamma^mu to the right, which gives us
    //
    //  -gamma^5 * adj(Dslash qOut) gamma^5 qIn
    //
    // Which is what we actually compute
    bilinear = g5 * adj(-Dslash_qOut) * g5 * qIn;
    compute_result(result.dslash_left);

    // Note: since gamma5^2 = 1 we can simplify
    bilinear = g5 * adj(-Dslash_qOut) * qIn;
    compute_result(result.dslash_gamma5_left);

    bilinear = g5 * adj(qOut) * g5 * Dslash_qIn;
    compute_result(result.dslash_right);

    bilinear = g5 * adj(qOut) * Dslash_qIn;
    compute_result(result.dslash_gamma5_right);

    bilinear = g5 * adj(qOut) * g5 * qIn;
    compute_result(result.scalar);

    bilinear = g5 * adj(qOut) * qIn;
    compute_result(result.psuedoscalar);

    LOG(Message) << "Complete. Writing results to " << par().output << std::endl;
    saveResult(par().output, "SubtractionOperators", result);
    envGet(HadronsSerializable, getName()) = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNPR_SubtractionOperators_hpp_
