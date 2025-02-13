/*
 * FermionFlow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MGradientFlow_FermionFlow_hpp_
#define Hadrons_MGradientFlow_FermionFlow_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Hadrons/Modules/MGradientFlow/Utils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Propagator Field Gradient Flow                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MGradientFlow)

class FermionFlowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FermionFlowPar,
                                    std::string, output,
                                    std::vector<std::string>, props,
                                    std::string, gauge,
                                    int, bc,
                                    int, steps,
                                    double, step_size,
                                    int, meas_interval,
                                    std::string, maxTau);
};

template <typename FImpl,typename GImpl,typename FlowAction>
class TFermionFlow: public Module<FermionFlowPar>
{
public:
    /*FERM_TYPE_ALIASES(FImpl,);
    INHERIT_GIMPL_TYPES(GImpl);*/
    BASIC_TYPE_ALIASES(FImpl,);
    GAUGE_TYPE_ALIASES(GImpl,);
    class GaugeResult : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GaugeResult,
                                        std::vector<double>,    plaquette,
                                        std::vector<double>,    rectangle,
                                        std::vector<double>,    clover,
                                        std::vector<double>,    topcharge,
                                        std::vector<double>,    action,
                                        std::vector<ComplexD>,  polyakov);
    };
public:
    // constructor
    TFermionFlow(const std::string name);
    // destructor
    virtual ~TFermionFlow(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // action
    FlowAction SG = FlowAction(3.0);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(WilsonFermionFlow,ARG(TFermionFlow<FIMPL,GIMPL,WilsonGaugeAction<GIMPL>>),MGradientFlow);

/******************************************************************************
 *                     TFermionFlow implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl,typename GImpl,typename FlowAction>
TFermionFlow<FImpl,GImpl,FlowAction>::TFermionFlow(const std::string name)
: Module<FermionFlowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl,typename GImpl,typename FlowAction>
std::vector<std::string> TFermionFlow<FImpl,GImpl,FlowAction>::getInput(void)
{
    std::vector<std::string> in = {par().gauge}; 
    for (std::string q : par().props) {
        in.push_back(q);
    }
    
    return in;
}

template <typename FImpl,typename GImpl,typename FlowAction>
std::vector<std::string> TFermionFlow<FImpl,GImpl,FlowAction>::getOutput(void)
{
    std::vector<std::string> out = {getName(),getName()+"_U"};
    for (int i = 1; i <= par().steps; i++) 
    {
        if ((i % par().meas_interval == 0) || (i == par().steps)) {
            std::stringstream st; st << i;
            for (int j = 0; j < par().props.size(); j++) {
                std::stringstream qt; qt << j;
                out.push_back(getName()+"_q"+qt.str()+"_"+st.str());
            }
        }
    }
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl,typename GImpl,typename FlowAction>
void TFermionFlow<FImpl,GImpl,FlowAction>::setup(void)
{
    envCreateLat(GaugeField, getName()+"_U");

    for (int j = 0; j < par().props.size(); j++) {
        std::stringstream qt; qt << j;
        envTmpLat(PropagatorField, "q"+qt.str()+"wf");
    }

    for (int i = 1; i <= par().steps; i++) 
    {
        if ((i % par().meas_interval == 0) || (i == par().steps)) {
            std::stringstream st; st << i;
            for (int j = 0; j < par().props.size(); j++) {
                std::stringstream qt; qt << j;
                envCreateLat(PropagatorField, getName()+"_q"+qt.str()+"_"+st.str());
            }
        }
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl,typename GImpl,typename FlowAction>
void TFermionFlow<FImpl,GImpl,FlowAction>::execute(void)
{
    std::string type = SG.action_name();
    std::string ga = "GaugeAction";
    std::string::size_type i = type.find(ga);
    if (i != std::string::npos) {
        type.erase(i, ga.length());
    }

    std::string props = "";
    for (std::string q : par().props) props += q + " ";
    LOG(Message) << "Setting up " << type << " Fermion Flow on '" << par().gauge << "' Gauge Field and "  
                 << props << ((par().props.size() > 1) ? "Fermion Propagators " : "Fermion Propagator ")
                 << "with ppp" << ((par().bc < 0) ? "a" : "p") << " boundary conditions and "
                 << par().steps << " step" << ((par().steps > 1) ? "s." : ".") << std::endl;

    std::vector<int> bc = {1,1,1};
    if (par().bc < 0) bc.push_back(-1);
    else bc.push_back(1);

    double mTau = -1.0;
    if(!par().maxTau.empty()) {
        LOG(Message) << "Using adaptive algorithm with maxTau = " << par().maxTau << std::endl;
        mTau = (double)std::stoi(par().maxTau);
    }

    auto &out     = envGet(HadronsSerializable, getName());
    auto &Uresult = out.template hold<GaugeResult>();

    Uresult.plaquette.resize(par().steps);
    Uresult.rectangle.resize(par().steps);
    Uresult.clover.resize(par().steps);
    Uresult.topcharge.resize(par().steps);
    Uresult.action.resize(par().steps);
    Uresult.polyakov.resize(par().steps);

    auto &U   = envGet(GaugeField, par().gauge);
    auto &Uwf = envGet(GaugeField, getName()+"_U");
    Uwf = U;

    for (int j = 0; j < par().props.size(); j++) {
        auto &qj = envGet(PropagatorField, par().props[j]);
        std::stringstream jt; jt << j;
        PropagatorField &qjwf = *env().template getObject<PropagatorField>(getName()+"_tmp_q"+jt.str()+"wf");
        qjwf = qj;
    }
    
    double time = 0;
    Evolution<FlowAction> evolve(3.0, par().step_size, mTau, par().step_size);
    if (mTau > 0) {
        unsigned int step = 0;
        do {
            step++;
            startTimer("gauge field step");
            std::vector<GaugeField> Wi = evolve.template evolve_gaugeFF_adaptive<GImpl,GaugeField,GaugeLinkField>(Uwf,bc);
            stopTimer("gauge field step");
            evolve.template gauge_status<GImpl,GaugeField,ComplexField,GaugeLinkField,GaugeResult>(Uwf,Uresult,step-1);
            if (step % par().meas_interval == 0) {
                std::stringstream st; st << step;
                for (int j = 0; j < par().props.size(); j++) {
                    std::stringstream jt; jt << j;
                    PropagatorField &qjwf = *env().template getObject<PropagatorField>(getName()+"_tmp_q"+jt.str()+"wf");
                    startTimer("fermion field "+jt.str()+" step");
                    evolve.template laplace_flow<PropagatorField,GImpl,GaugeField,GaugeLinkField>(Wi[0],Wi[1],Wi[2],qjwf);
                    stopTimer("fermion field "+jt.str()+" step");
                    auto &qji = envGet(PropagatorField, getName()+"_q"+jt.str()+"_"+st.str());
                    qji = qjwf;
                }
            }
        } while (evolve.taus < mTau);
    } else {
        for (unsigned int step = 1; step <= par().steps; step++) {
            startTimer("gauge field step");
            std::vector<GaugeField> Wi = evolve.template evolve_gaugeFF<GImpl,GaugeField,GaugeLinkField>(Uwf,bc);
            stopTimer("gauge field step");
            evolve.template gauge_status<GImpl,GaugeField,ComplexField,GaugeLinkField,GaugeResult>(Uwf,Uresult,step-1);
            if ((step % par().meas_interval == 0) || (step == par().steps)) {
                std::stringstream st; st << step;
                for (int j = 0; j < par().props.size(); j++) {
                    std::stringstream jt; jt << j;
                    PropagatorField &qjwf = *env().template getObject<PropagatorField>(getName()+"_tmp_q"+jt.str()+"wf");
                    startTimer("fermion field "+jt.str()+" step");
                    evolve.template laplace_flow<PropagatorField,GImpl,GaugeField,GaugeLinkField>(Wi[0],Wi[1],Wi[2],qjwf);
                    stopTimer("fermion field "+jt.str()+" step");
                    auto &qji = envGet(PropagatorField, getName()+"_q"+jt.str()+"_"+st.str());
                    qji = qjwf;
                }
            }
        }
    }
    saveResult(par().output,"gauge_obs",Uresult);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGradientFlow_FermionFlow_hpp_
