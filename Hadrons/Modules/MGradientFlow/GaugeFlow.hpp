/*
 * GaugeFlow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGradientFlow_GaugeFlow_hpp_
#define Hadrons_MGradientFlow_GaugeFlow_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MGradientFlow/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                        Gauge Field Gradient Flow                           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGradientFlow)

class GaugeFlowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeFlowPar,
                                    std::string, output,
                                    std::string, gauge,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval,
                                    std::string, maxTau,
                                    std::string, c_plaq,
                                    std::string, c_rect); 
};

template <typename GImpl,typename FlowAction>
class TGaugeFlow: public Module<GaugeFlowPar>
{
public:
    INHERIT_GIMPL_TYPES(GImpl);
    class Result : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        std::vector<double>,    plaquette,
                                        std::vector<double>,    rectangle,
                                        std::vector<double>,    clover,
                                        std::vector<double>,    topcharge,
                                        std::vector<double>,    action,
                                        std::vector<ComplexD>,  polyakov);
    };
public:
    // constructor
    TGaugeFlow(const std::string name);
    // destructor
    virtual ~TGaugeFlow(void) {};
    // run flow
    virtual void runGaugeFlow(Evolution<FlowAction> &evolve, double &mTau);
    virtual void setupGaugeFlow(FlowAction &SG);
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonFlow, ARG(TGaugeFlow<GIMPL,WilsonGaugeAction<GIMPL>>), MGradientFlow);
MODULE_REGISTER_TMP(SymanzikFlow, ARG(TGaugeFlow<GIMPL,SymanzikGaugeAction<GIMPL>>), MGradientFlow);
MODULE_REGISTER_TMP(ZeuthenFlow, ARG(TGaugeFlow<GIMPL,ZeuthenGaugeAction<GIMPL>>), MGradientFlow);
MODULE_REGISTER_TMP(CustomFlow, ARG(TGaugeFlow<GIMPL,PlaqPlusRectangleAction<GIMPL>>), MGradientFlow);

/******************************************************************************
 *                     TGaugeFlow implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
TGaugeFlow<GImpl,FlowAction>::TGaugeFlow(const std::string name)
: Module<GaugeFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
std::vector<std::string> TGaugeFlow<GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};
    
    return in;
}

template <typename GImpl,typename FlowAction>
std::vector<std::string> TGaugeFlow<GImpl,FlowAction>::getOutput(void)
{
    std::vector<std::string> out = {getName(),getName()+"_U"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::setup(void)
{
    envCreateLat(GaugeField, getName()+"_U");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// run flow ////////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::runGaugeFlow(Evolution<FlowAction> &evolve, double &mTau)
{
    auto &out    = envGet(HadronsSerializable, getName());
    auto &result = out.template hold<Result>();

    auto &U   = envGet(GaugeField, par().gauge);
    auto &Uwf = envGet(GaugeField, getName()+"_U");
    Uwf = U;

    if (par().steps == 0) { 
        // if steps = 0, give the status of gauge field without flowing
        result.plaquette.resize(1);
        result.rectangle.resize(1);
        result.clover.resize(1);
        result.topcharge.resize(1);
        result.action.resize(1);
        result.polyakov.resize(1);
        evolve.template gauge_status<GImpl,GaugeField,ComplexField,GaugeLinkField,Result>(Uwf,result,0); 
    } else {
        result.plaquette.resize(par().steps);
        result.rectangle.resize(par().steps);
        result.clover.resize(par().steps);
        result.topcharge.resize(par().steps);
        result.action.resize(par().steps);
        result.polyakov.resize(par().steps);
        if (mTau > 0) {
            unsigned int step = 0;
            do {
                step++;
                evolve.template evolve_gauge_adaptive<GImpl,GaugeField>(Uwf);
                if (step % par().meas_interval == 0) {
                    evolve.template gauge_status<GImpl,GaugeField,ComplexField,GaugeLinkField,Result>(Uwf,result,step-1);
                }
            } while (evolve.taus < mTau);
        } else {
            for (unsigned int step = 1; step <= par().steps; step++) {
                evolve.template evolve_gauge<GImpl,GaugeField>(Uwf);
                if (step % par().meas_interval == 0) {
                    evolve.template gauge_status<GImpl,GaugeField,ComplexField,GaugeLinkField,Result>(Uwf,result,step-1);
                }
            }
        }
    }
    saveResult(par().output,"gauge_obs",result);
}

template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::setupGaugeFlow(FlowAction &SG)
{
    std::string type = SG.action_name();
    LOG(Message) << "Setting up " << type << " Flow on '" << par().gauge << "' with " << par().steps
                 << " step" << ((par().steps != 1) ? "s." : ".") << std::endl;

    double mTau = -1.0;
    if(!par().maxTau.empty()) {
        LOG(Message) << "Using adaptive algorithm with maxTau = " << par().maxTau << std::endl;
        mTau = std::stod(par().maxTau);
    }

    if constexpr (std::is_same_v<FlowAction, PlaqPlusRectangleAction<GImpl>>) {
        Evolution<FlowAction> evolve(std::stod(par().c_plaq), std::stod(par().c_rect), par().step_size, mTau, par().step_size);
        runGaugeFlow(evolve,mTau);
    } else {
        Evolution<FlowAction> evolve(3.0, par().step_size, mTau, par().step_size);
        runGaugeFlow(evolve,mTau);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl,typename FlowAction>
void TGaugeFlow<GImpl,FlowAction>::execute(void)
{
    // create action -> if c_plaq and c_rect are used, called PlaqPlusRectangleAction
    if constexpr (std::is_same_v<FlowAction, PlaqPlusRectangleAction<GImpl>>) {
        if (par().c_plaq.empty() || par().c_plaq.empty()) {
            std::cerr << "Error: to use PlaqPlusRectangleAction (CustomFlow), pass some value to both c_plaq and c_rect." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        FlowAction SG = FlowAction(std::stod(par().c_plaq),std::stod(par().c_rect));
        LOG(Message) << SG.LogParameters();
        setupGaugeFlow(SG);
    } else {
        FlowAction SG = FlowAction(3.0);
        setupGaugeFlow(SG);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_GaugeFlow_hpp_
