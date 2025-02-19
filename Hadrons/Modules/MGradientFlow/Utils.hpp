/*
 * Utils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2022
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Matthew Black    <matthewkblack@protonmail.com>
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
#ifndef Hadrons_MGradientFlow_Utils_hpp_
#define Hadrons_MGradientFlow_Utils_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MGradientFlow)

// additional action(s) /////////////////////////////////////////////////////
template <class GImpl>
class ZeuthenGaugeAction {
public:
    INHERIT_GIMPL_TYPES(GImpl);
    
    double beta;
    SymanzikGaugeAction<GImpl> SG;

    ZeuthenGaugeAction(double b): beta(b),SG(SymanzikGaugeAction<GImpl>(b)) {};

    virtual std::string action_name(){return "ZeuthenGaugeAction";}

    virtual double S(const GaugeField &U) {
        return SG.S(U);
    };

    virtual void deriv(const GaugeField &Umu, GaugeField &dSdU) {
                                                //  beta = 3.0, cl = -1.0/12.0 -> Symanzik
        double factor_p = 5.0/double(Nc)*0.5;   //   5.0 = beta*(1.0-8.0*cl)
        double factor_r = -0.25/double(Nc)*0.5; // -0.25 = beta*cl

        GridBase *grid = Umu.Grid();

        std::vector<GaugeLinkField> U (Nd,grid);
        std::vector<GaugeLinkField> U2(Nd,grid);

        for(int mu=0;mu<Nd;mu++){
            U[mu] = PeekIndex<LorentzIndex>(Umu,mu);
            WilsonLoops<GImpl>::RectStapleDouble(U2[mu],U[mu],mu);
        }

        GaugeLinkField dSdU_mu(grid);
        GaugeLinkField staple(grid);
        GaugeLinkField tmp(grid),tmq(grid),tmr(grid);

        for (int mu=0; mu < Nd; mu++){
            // Staple in direction mu
            WilsonLoops<GImpl>::Staple(staple,Umu,mu);
            tmp = Ta(U[mu]*staple)*factor_p;

            WilsonLoops<GImpl>::RectStaple(Umu,staple,U2,U,mu);
            tmp = tmp + Ta(U[mu]*staple)*factor_r;

            tmq = (adj(Cshift(U[mu],mu,-1)) * Cshift(tmp,mu,-1) * Cshift(U[mu],mu,-1));
            tmr = (U[mu] * Cshift(tmp,mu,1) * adj(U[mu]));

            dSdU_mu = 5.0/6.0*tmp + 1.0/12.0*tmq + 1.0/12.0*tmr;
            PokeIndex<LorentzIndex>(dSdU, dSdU_mu, mu);
        }
    };
};

// clover //////////////////////////////////////////////////////////////////////
template <typename GImpl, typename ComplexField, typename GaugeLorentz, typename GaugeMat>
void siteClover(ComplexField &Clov, const GaugeLorentz &U)
{
    GaugeMat Fmn(U.Grid()), Cmn(U.Grid()), scaledUnit(U.Grid()), Umu(U.Grid());
    Clov = Zero();
    for (int mu = 1; mu < Nd; mu++) {
        for (int nu = 0; nu < mu; nu++) {
            Umu = PeekIndex<LorentzIndex>(U, mu);
            scaledUnit = (1.0/Nc) * (adj(Umu) * Umu);
            WilsonLoops<GImpl>::FieldStrength(Fmn, U, mu, nu);
            Cmn = Fmn - trace(Fmn) * scaledUnit;
            Clov = Clov - trace(Cmn * Cmn);
        }
    }
}

template <typename GImpl, typename ComplexField, typename GaugeLorentz, typename GaugeMat>
double avgClover(const GaugeLorentz &Umu) 
{
    ComplexField Clov(Umu.Grid());

    siteClover<GImpl,ComplexField,GaugeLorentz,GaugeMat>(Clov, Umu);
    auto Tc = sum(Clov);
    auto c = TensorRemove(Tc);

    double vol = Umu.Grid()->gSites();

    return c.real() / vol;
}

// field evolution /////////////////////////////////////////////////////////////
template <typename FlowAction>
class Evolution {
    public:
        double epsilon, maxTau, taus;
        FlowAction SG;
        Evolution(double beta, double step, double mTau, double ts) : 
            SG(FlowAction(beta)), epsilon(step), maxTau(mTau), taus(ts) {};

        template <typename GImpl,typename GaugeField>
        std::vector<GaugeField> gauge_RK(GaugeField U) {

            std::vector<GaugeField> Wi;

            GaugeField Z(U.Grid());
            GaugeField tmp(U.Grid());
            Wi.push_back(U);                            // W0
            SG.deriv(U, Z);                                
            Z *= 0.25;                                  // Z0 = 1/4 * F(U)
            GImpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0
            Wi.push_back(U);                            // W1

            Z *= -17.0/8.0;
            SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 + Z1
            Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
            GImpl::update_field(Z, U, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1
            Wi.push_back(U);                            // W2

            Z *= -4.0/3.0;
            SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) + Z2
            Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
            GImpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2
            Wi.push_back(U);                            // W3

            return Wi;
        };

        template <typename GImpl,typename GaugeField>
        std::vector<GaugeField> gauge_RK_adaptive(GaugeField U) {

            std::vector<GaugeField> Wi;

            if (maxTau - taus < epsilon){
                epsilon = maxTau-taus;
            }
            GaugeField Z(U.Grid());
            GaugeField Zprime(U.Grid());
            GaugeField tmp(U.Grid()), Uprime(U.Grid());
            Uprime = U;
            Wi.push_back(U);                            // W0
            SG.deriv(U, Z);
            Zprime = -Z;
            Z *= 0.25;                                  // Z0 = 1/4 * F(U)
            GImpl::update_field(Z, U, -2.0*epsilon);    // U = W1 = exp(ep*Z0)*W0
            Wi.push_back(U);                            // W1

            Z *= -17.0/8.0;
            SG.deriv(U, tmp); Z += tmp;                 // -17/32*Z0 +Z1
            Zprime += 2.0*tmp;
            Z *= 8.0/9.0;                               // Z = -17/36*Z0 +8/9*Z1
            GImpl::update_field(Z, U, -2.0*epsilon);    // U_= W2 = exp(ep*Z)*W1
            Wi.push_back(U);                            // W2

            Z *= -4.0/3.0;
            SG.deriv(U, tmp); Z += tmp;                 // 4/3*(17/36*Z0 -8/9*Z1) +Z2
            Z *= 3.0/4.0;                               // Z = 17/36*Z0 -8/9*Z1 +3/4*Z2
            GImpl::update_field(Z, U, -2.0*epsilon);    // V(t+e) = exp(ep*Z)*W2      
            Wi.push_back(U);                            // W3
            
            GImpl::update_field(Zprime, Uprime, -2.0*epsilon); // V'(t+e) = exp(ep*Z')*W0
            Wi.push_back(Uprime);                       // Uprime

            return Wi;
        };

        template <typename GaugeField>
        void adaptive_eps(const GaugeField& U, const GaugeField& Uprime) {
            // Compute distance as norm^2 of the difference
            GaugeField diffU = U - Uprime;
            double diff = norm2(diffU);   
            // adjust integration step  

            taus += epsilon;
            epsilon = epsilon*0.95*std::pow(1e-4/diff,1./3.);
        };

        template <typename GImpl,typename GaugeField>
        void evolve_gauge(GaugeField &U) {
            std::vector<GaugeField> Wi = gauge_RK<GImpl,GaugeField>(U);
            U = Wi[3];
        };

        template <typename GImpl,typename GaugeField>
        void evolve_gauge_adaptive(GaugeField &U) {
            std::vector<GaugeField> Wi = gauge_RK_adaptive<GImpl,GaugeField>(U);
            adaptive_eps(Wi[3],Wi[4]);
            U = Wi[3];
        };

        template <typename GaugeField,typename GaugeLinkField>
        void gauge_apply_boundary(GaugeField &Umu, std::vector<int> bc) {
            GaugeLinkField tmp1(Umu.Grid());
            GaugeLinkField tmp2(Umu.Grid());
            GaugeLinkField tmp3(Umu.Grid());
            Lattice<iScalar<vInteger>> coord(Umu.Grid());

            for (int mu = 0; mu < Nd; mu++) {
                LatticeCoordinate(coord,mu);
                
                tmp1 = PeekIndex<LorentzIndex>(Umu,mu);
                tmp2 = (double)bc[mu]*tmp1;
                int dimSize = Umu.Grid()->GlobalDimensions()[mu] - 1;
                tmp3 = where((coord == dimSize), tmp2, tmp1);
                PokeIndex<LorentzIndex>(Umu, tmp3, mu);
            }
        };

        template <typename FImpl,typename GImpl,typename GaugeField,typename GaugeLinkField>
        FImpl generic_laplace(double a, double b, GaugeField &Umu, const FImpl& x_in, int skip_axis) {
            double Nx = Nd;
            if (skip_axis != -1) Nx--;

            FImpl x_out = (a + -2.0*Nx*b) * x_in;
            for (int mu = 0; mu < Nd; mu++) {
                if (mu != skip_axis) {
                    GaugeLinkField U = PeekIndex<LorentzIndex>(Umu, mu);
                    x_out += b*(GImpl::CovShiftForward(U,mu,x_in) + GImpl::CovShiftBackward(U,mu,x_in));
                }
            }
            return x_out;
        };

        template <typename FImpl,typename GImpl,typename GaugeField,typename GaugeLinkField>
        void laplace_flow(GaugeField &W0, GaugeField &W1, GaugeField &W2, FImpl &prop) {
            FImpl psi1 = prop + (epsilon/4.0)*generic_laplace<FImpl,GImpl,GaugeField,GaugeLinkField>(0.0, 1.0, W0, prop, -1);
            FImpl psi2 = prop + (8.0*epsilon/9.0)*generic_laplace<FImpl,GImpl,GaugeField,GaugeLinkField>(0.0, 1.0, W1, psi1, -1) - (2.0*epsilon/9.0)*generic_laplace<FImpl,GImpl,GaugeField,GaugeLinkField>(0.0, 1.0, W0, prop, -1);
            FImpl psi3 = psi1 + (3.0*epsilon/4.0)*generic_laplace<FImpl,GImpl,GaugeField,GaugeLinkField>(0.0, 1.0, W2, psi2, -1);

            prop = psi3;
        };

        template <typename GImpl,typename GaugeField,typename GaugeLinkField>
        std::vector<GaugeField> evolve_gaugeFF(GaugeField &U, std::vector<int> &bc) {
            std::vector<GaugeField> Wi = gauge_RK<GImpl,GaugeField>(U);
            U = 1.0*Wi[3];

            gauge_apply_boundary<GaugeField,GaugeLinkField>(Wi[0],bc);
            gauge_apply_boundary<GaugeField,GaugeLinkField>(Wi[1],bc);
            gauge_apply_boundary<GaugeField,GaugeLinkField>(Wi[2],bc);

            return Wi;
        };

        // gauge field status //////////////////////////////////////////////////////////
        template <typename GImpl,typename GaugeField,typename ComplexField,typename GaugeLinkField,typename Result>
        void gauge_status(GaugeField &Umu, Result &result, int index)
        {
            double Q = WilsonLoops<GImpl>::TopologicalCharge(Umu);
            double plaq = WilsonLoops<GImpl>::avgPlaquette(Umu);
            double rect = WilsonLoops<GImpl>::avgRectangle(Umu);
            double clov = avgClover<GImpl,ComplexField,GaugeField,GaugeLinkField>(Umu);
            double act = SG.S(Umu);
            ComplexD poly = WilsonLoops<GImpl>::avgPolyakovLoop(Umu);

            result.plaquette[index] = plaq;
            result.rectangle[index] = rect;
            result.clover[index]    = clov;
            result.topcharge[index] =    Q;
            result.action[index]    =  act;
            result.polyakov[index]  = poly;
        };
};

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_Utils_hpp_
