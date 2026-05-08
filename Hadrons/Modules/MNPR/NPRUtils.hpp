/*
 * NPRUtils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
 * Author: Felix Erben <felix.erben@ed.ac.uk>
 * Author: Felix Erben <ferben@ed.ac.uk>
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

#ifndef Hadrons_MNPR_NPRUtils_hpp_
#define Hadrons_MNPR_NPRUtils_hpp_

#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MNPR)

template <typename FImpl>
class NPRUtils
{
public:
    FERM_TYPE_ALIASES(FImpl,)
private:
    //this is static so that we can use the cached momentum space windows without recomputing for each config
    static std::map<std::string, ComplexField>& get_cache(){
	static std::map<std::string, ComplexField> cache; 
	return cache; 
	};
    //common routine for the FFT and caching for any window function
    static ComplexField& cache_and_return(std::string key, ComplexField& w_x);
public:
    //v1
    static SpinColourSpinColourMatrix tensorProdSum(PropagatorField &tsum, PropagatorField &a, PropagatorField &b);
    //v2
    static SpinColourSpinColourMatrix tensorProdSum(LatticeColourMatrix &tmp_colour, PropagatorField &tsum, PropagatorField &a, PropagatorField &b);
    //v3
    static SpinColourSpinColourMatrix tensorProdSum_v3(PropagatorField &tsum, PropagatorField &a, PropagatorField &b);

    static void tensorSiteProd(SpinColourSpinColourMatrix &lret, SpinColourMatrixScalar &a, SpinColourMatrixScalar &b);        
    static void tensorSiteProd(vSpinColourSpinColourMatrix &lret, vSpinColourMatrix &a, vSpinColourMatrix &b);
    static void tensorProd(LatticeSpinColourSpinColourMatrix &res, PropagatorField &a, PropagatorField &b); 
    // covariant derivative
    static void dslash(PropagatorField &in, const PropagatorField &out,
        const GaugeField &Umu);
    static void phase(ComplexField &bilinearPhase, std::vector<Real> pIn, std::vector<Real> pOut);
    static void dot(ComplexField &pDotX, std::vector<Real> p);
    static void getSymmetricCoords(std::vector<ComplexField> &sym_coords, GridBase *grid);
    static const ComplexField& getRectWindow(GridBase *grid); //trivial window for debugging
    static const ComplexField& getGaussianWindow(GridBase *grid, double sigma);
};


// Tensor product of two PropagatorFields (Lattice Spin Colour Matrices in many FImpls)
// v1
template <typename FImpl>
SpinColourSpinColourMatrix NPRUtils<FImpl>::tensorProdSum(PropagatorField &tsum, PropagatorField &a, PropagatorField &b)
{
    GridBase *grid = a.Grid();
    SpinColourSpinColourMatrix result;
    double t_acc=0.0;
    double t_sum=0.0;

    for(int si=0; si < Ns; ++si)
	{
        for(int sj=0; sj < Ns; ++sj)
	    {
            for (int ci=0; ci < Nc; ++ci)
	        {
                for (int cj=0; cj < Nc; ++cj)
	            {
                    t_acc -= usecond();    
                    tsum = peekColour(peekSpin(a, si, sj), ci, cj) * b;
                    t_acc += usecond();    
		    t_sum -= usecond();
                    result()(si,sj)(ci,cj) = sum_large(tsum)();
		    t_sum += usecond();
                }
            }
        }
    }
    
    uint64_t bytes  = Nc*Nc*Ns*Ns*3*sizeof(vSpinColourMatrix)*grid->oSites();
    uint64_t sbytes = Nc*Nc*Ns*Ns*sizeof(vSpinColourMatrix)*grid->oSites();
    std::cout << " tacc "<< t_acc << " tsum " << t_sum << " us "<<std::endl;
    std::cout << " bytes "<< bytes <<std::endl;
    std::cout << "acc "<< (double)bytes/t_acc << " MB/s "<<std::endl;
    std::cout << "sum "<< (double)sbytes/t_sum << " MB/s "<<std::endl;
    return result;
}

//v2. just adds an additional temporary such that we don't waste spin index iterations
template <typename FImpl>
SpinColourSpinColourMatrix NPRUtils<FImpl>::tensorProdSum(LatticeColourMatrix &tmp_colour, PropagatorField &tsum, PropagatorField &a, PropagatorField &b)
{
    GridBase *grid = a.Grid(); 

    SpinColourSpinColourMatrix result;
    double t_acc=0.0;
    double t_sum=0.0;

    for(int si=0; si < Ns; ++si)
	{
        for(int sj=0; sj < Ns; ++sj)
	    {
            tmp_colour = peekSpin(a, si, sj); 
	    
 	    for (int ci=0; ci < Nc; ++ci)
	        {
                for (int cj=0; cj < Nc; ++cj)
	            {
                    t_acc -= usecond();    
		    tsum = peekColour(tmp_colour, ci, cj) * b;
                    t_acc+=usecond();
		    t_sum-=usecond();
		    result()(si,sj)(ci,cj) = sum_large(tsum)();
                    t_sum+=usecond();
		}
            }
        }
    } 
    uint64_t bytes  = Nc*Nc*Ns*Ns*3*sizeof(vSpinColourMatrix)*grid->oSites();
    uint64_t sbytes = Nc*Nc*Ns*Ns*sizeof(vSpinColourMatrix)*grid->oSites();
    std::cout << " tacc "<< t_acc << " tsum " << t_sum << " us "<<std::endl;
    std::cout << " bytes "<< bytes <<std::endl;
    std::cout << "acc "<< (double)bytes/t_acc << " MB/s "<<std::endl;
    std::cout << "sum "<< (double)sbytes/t_sum << " MB/s "<<std::endl;
    return result;
}

//v3.
template <typename FImpl>
SpinColourSpinColourMatrix NPRUtils<FImpl>::tensorProdSum_v3(PropagatorField &tsum,PropagatorField &a, PropagatorField &b)
{
    vComplex dummy;
    SpinColourSpinColourMatrix result;
    GridBase *grid = a.Grid();
    int Nsimd = grid->Nsimd();
    double t_acc=0.0;
    double t_sum=0.0;
    autoView(tsum_v,tsum,AcceleratorWrite);
    autoView(a_v,a,AcceleratorRead);
    autoView(b_v,b,AcceleratorRead);
    for(int si=0; si < Ns; ++si){
    for(int sj=0; sj < Ns; ++sj){
    for(int ci=0; ci < Nc; ++ci){ 
    for(int cj=0; cj < Nc; ++cj){
        t_acc -= usecond();
        accelerator_for(ss,grid->oSites(),Nsimd, {
            typedef decltype(coalescedRead(vTComplex())) calcComplex;
            calcComplex aa;
            aa()()() = a_v(ss)()(si,sj)(ci,cj);
            auto bb = b_v(ss);
            coalescedWrite(tsum_v[ss],aa*bb);        
        });
        t_acc += usecond();
        t_sum -= usecond();
        result()(si,sj)(ci,cj) = sum_large(tsum)();
        t_sum += usecond();
    }}}}
    uint64_t bytes  = Ns*Ns*Nc*Nc*3*sizeof(vSpinColourMatrix)*grid->oSites();
    uint64_t sbytes = Ns*Ns*Nc*Nc*sizeof(vSpinColourMatrix)*grid->oSites();
    std::cout << " tacc "<< t_acc << " tsum " << t_sum << " us "<<std::endl;
    std::cout << " bytes "<< bytes <<std::endl;
    std::cout << "acc "<< (double)bytes/t_acc << " MB/s "<<std::endl;
    std::cout << "sum "<< (double)sbytes/t_sum << " MB/s "<<std::endl;
    return result;
}

    
// Tensosumroduct on a ssingle sisumonly
template <typename FImpl>
void NPRUtils<FImpl>::tensorSiteProd(SpinColourSpinColourMatrix &lret, SpinColourMatrixScalar &a, SpinColourMatrixScalar &b)
{
    for(int si=0; si < Ns; ++si)
    {
    for(int sj=0; sj < Ns; ++sj)
    {
        for (int ci=0; ci < Nc; ++ci)
	{
        for (int cj=0; cj < Nc; ++cj)
	{
            const ComplexD val = TensorRemove(a()(si,sj)(ci,cj));
            lret()(si,sj)(ci,cj) = val * b();
        }}
    }}
}


// Computes gamma^mu D_mu for the given input field. Currently uses the
// symmetric derivative, though this could change in the future.
template <typename FImpl>
void NPRUtils<FImpl>::dslash(PropagatorField &out, const PropagatorField &in,
        const GaugeField &Umu)
{
    assert(&out != &in);
    out = Zero();
    PropagatorField tmp(Umu.Grid());
    typename FImpl::GaugeLinkField U(Umu.Grid());
    for (int mu = 0; mu < Nd; mu++)
    {
        // Overall formula:
        // tmp(x) = U_\mu(x) in(x + \hat{\mu}) - U_\mu^\dag(x - \hat{\mu}) in(x - \hat{\mu})
        U = peekLorentz(Umu, mu);
        tmp = FImpl::CovShiftForward(U, mu, in);
        tmp = tmp - FImpl::CovShiftBackward(U, mu, in);

        Gamma gamma_mu = Gamma::gmu[mu];
        out += gamma_mu * tmp;
    }
    out = 0.5 * out;
}


//// Compute phases for phasing propagators
// bilinearPhase = exp(-i (pIn - pOut) \cdot x)
template <typename FImpl>
void NPRUtils<FImpl>::phase(ComplexField &bilinearPhase, std::vector<Real> pIn, std::vector<Real> pOut)
{
    bilinearPhase = Zero();
    ComplexField coordinate(bilinearPhase.Grid());
    Coordinate                  latt_size = GridDefaultLatt();
    for (int mu = 0; mu < Nd; mu++)
    {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;

        bilinearPhase += coordinate * (pIn[mu] - pOut[mu]);
    }
    Complex Ci = Complex(0.0, 1.0);
    bilinearPhase = exp(-Ci * bilinearPhase);
}


// pDotX = p \cdot x
template <typename FImpl>
void NPRUtils<FImpl>::dot(ComplexField &pDotX, std::vector<Real> p)
{
    ComplexField coordinate(pDotX.Grid());
    Coordinate                  latt_size = GridDefaultLatt();
    pDotX = Zero();
    for (int mu = 0; mu < Nd; mu++)
    {
        LatticeCoordinate(coordinate, mu);
        coordinate = (2 * M_PI / latt_size[mu]) * coordinate;
        pDotX += coordinate * p[mu];
    }
}

/*
class WindowCache {
    private:
	    //this is static so that we can use the cached momentum space windows without recomputing for each config
	    static std::map<std::string, ComplexField>& get_cache(){
	    	static std::map<std::string, ComplexField> cache; 
		return cache; 
	    };
	    //common routine for the FFT and caching for any window function
	    static ComplexField& cache_and_return(std::string key, ComplexField& w_x);
    public:
	    //various window functions. 
	    static const ComplexField& getGaussianWindow(ComplexField& coordinate, double sigma,  double N);
	    static const ComplexField& getCosineWindow(ComplexField& coordinate, double a0, double a1, double a2, double alpha, double N);
	    static const ComplexField& getExpWindow(ComplexField& coordinate, double tau, double N);

};

*/

template <typename FImpl>
typename NPRUtils<FImpl>::ComplexField &
NPRUtils<FImpl>::cache_and_return(std::string key, ComplexField& w_x)
{
    auto& cache = get_cache();
    GridBase *grid = w_x.Grid();

    FFT fft((GridCartesian *)grid);

    //emplace returns a pair, the first being an interator to the element. This gets around the fact that ComplexField has no default constructor
    //which the map does not like when it looks for an element that hasn't been created yet with a particular key. 
    auto insertion = cache.emplace(std::piecewise_construct, std::forward_as_tuple(key),std::forward_as_tuple(grid));

    ComplexField &w_p = insertion.first->second; 

    fft.FFT_all_dim(w_p, w_x, FFT::forward);

    //momentum space site sum
    LOG(Message) << "momentum space window site sum: " << sum(w_p) << std::endl;

    return w_p; 

}

template <typename FImpl>
void NPRUtils<FImpl>::getSymmetricCoords(std::vector<ComplexField> &sym_coords, GridBase *grid)
{

    //cache the symmetric coordinate vector to use for the whole run, we don't want to compute it every time. 
    static std::vector<ComplexField> coord_cache;

    //if the cache is not populated, we populate it
    if (coord_cache.empty()) {

	int ndim = grid->_ndimension;
	coord_cache.resize(ndim, ComplexField(grid));

	LatticeInteger dist(grid); 
	Coordinate latt_size = GridDefaultLatt();


	for(int mu=0;mu<grid->_ndimension;mu++){
    	    LatticeCoordinate(dist, mu);
            LatticeCoordinate(coord_cache[mu],mu);
	
	    const int L = latt_size[mu];
	    const int Lhalf = L/2;

	    coord_cache[mu] = where((dist<Integer(Lhalf)),coord_cache[mu],coord_cache[mu]-Real(L));	

    	}

    }

    sym_coords = coord_cache; 

}

//rectangular window, we set every site to 1 in position space. This is effectively no window, but I am using it for testing purposes. 
template <typename FImpl>
const typename NPRUtils<FImpl>::ComplexField &
NPRUtils<FImpl>::getRectWindow(GridBase *grid)
{
    auto& cache = get_cache();
    std::string key = "rectangle";

    auto it = cache.find(key);
    if (it != cache.end()){
        return it->second;
    }

    ComplexField w_x(grid);
    w_x = 1.0; //this sets every site to 1.0

    //checksum, should be lattice volume in lattice units 
    LOG(Message) <<"Window site sum: " << sum(w_x) << std::endl;

    return cache_and_return(key, w_x);
}

template <typename FImpl>
const typename NPRUtils<FImpl>::ComplexField & 
NPRUtils<FImpl>::getGaussianWindow(GridBase* grid, double fwhm)
{
    auto& cache = get_cache();
    std::string key = "gaussian_" + std::to_string(fwhm);
	
    //iterator to the element labelled by 'key'
    auto it = cache.find(key);
    if (it != cache.end()){
	return it->second; //return the cache entry if it exists
    }

    //we only build the window function if it doesn't exist in the cache already. 
    
    ComplexField w_x(grid);
    ComplexField rsq(grid); //holds x \cdot x
    rsq = Zero();

    const int dim = grid->_ndimension; 
    Coordinate latt_size = GridDefaultLatt(); 
    std::vector<double> fact(dim); //holds the constant factor multiplying x^2  	
    std::vector<ComplexField> xmu;
    //xmu.resize(dim, RealField(grid));

    int l_x = latt_size[0];
    double sigma = fwhm / (2.0 * sqrt(2.0 * log(2.0))*l_x); 

    getSymmetricCoords(xmu, grid);

    for(int mu=0; mu<dim; mu++){
	fact[mu] = (-0.5/((sigma*sigma)*(l_x*l_x)));
	rsq += xmu[mu] * xmu[mu] * fact[mu]; 
    }

    w_x = exp(rsq);
    
    //normalization
    //Real volInv = 1.0/TensorRemove(sum(w_x)).real();
    //w_x = w_x * volInv;
	
    //output the sum over all sites of the position space window
    LOG(Message) <<"Window site sum: " << sum(w_x) << std::endl;

    //check window value at zero
    Complex val0;
    std::vector<int> z(dim, 0);
    peekSite(val0, w_x, z);

    LOG(Message) << "w_x(0) = " << val0 << std::endl; 


    return cache_and_return(key, w_x); 
    
}


//passing in a string for the window in the xml, this is an efficient way to classify 
//which window we want to use, and assign local variables to such a window at runtime.
enum class WindowType{NONE, RECTANGLE, GAUSSIAN, DELTA, ZERO, STEP}; 

inline WindowType parseWindowType(const std::string &s){
    static const std::unordered_map<std::string, WindowType> windowMap{
        {"none", WindowType::NONE},
	{"rectangle", WindowType::RECTANGLE},
	{"gaussian", WindowType::GAUSSIAN},
	{"delta", WindowType::DELTA},
	{"zero", WindowType::ZERO},
	{"step", WindowType::STEP}
    };

    auto it = windowMap.find(s);
    if (it == windowMap.end()){
	LOG(Warning) << "Unknown window type " << s << "defaulting to no window, type NONE" << std::endl;
        return WindowType::NONE;
    }
    
    return it->second; 
    
}


END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif
