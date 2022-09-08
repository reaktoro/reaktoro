// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "EquilibriumSpecs.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace {

/// Return a Species object in a Database with given formula and aggregate state
auto getSpecies(Database const& db, String const& formula, AggregateState aggstate) -> Species
{
    const auto selected = db.speciesWithAggregateState(aggstate);
    const auto idx = selected.findWithFormula(formula);
    if(idx < selected.size()) return selected[idx];
    else return Species();
}

/// Return a Species object in a Database with given formula and aqueous aggregate state.
auto getAqueousSpecies(Database const& db, String const& formula) -> Species
{
    return getSpecies(db, formula, AggregateState::Aqueous);
}

/// Return a Species object in a Database with given formula and gaseous aggregate state.
auto getGaseousSpecies(Database const& db, String const& formula) -> Species
{
    return getSpecies(db, formula, AggregateState::Gas);
}

/// Assemble the complete system of equation constraints to be satisfied at chemical equilibrium.
auto assembleEquationConstraints(EquilibriumSpecs const& specs) -> EquationConstraints
{
    // Create copies of all provided equation constraints to avoid full capture of `specs` in lambda function below
    const auto econstraints_single = specs.equationConstraintsSingle();
    const auto econstraints_system = specs.equationConstraintsSystem();

    // The complete system of equation constraints that will be created below.
    EquationConstraints econstraints;

    // Collect the id constraints across the single equation constraints
    for(auto const& x : econstraints_single)
        econstraints.ids.push_back(x.id);

    // Collect the id constraints across the system of equation constraints
    for(auto const& x : econstraints_system)
        econstraints.ids.insert(econstraints.ids.end(), x.ids.begin(), x.ids.end());

    // The total number of equation constraints
    const auto num_econstraints = econstraints.ids.size();

    // Create the final constraint vector function
    econstraints.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w) -> VectorXr
    {
        // Define auxiliary offset variable to keep track of the entries in `res` to be filled in
        auto offset = 0;

        // The vector with the resulting evaluation of all equation constraints
        VectorXr res(num_econstraints);

        // Evaluate all single equation constraints
        for(auto const& x : econstraints_single)
            res[offset++] = x.fn(state, p, w);

        // Evaluate all system of equation constraints
        for(auto const& x : econstraints_system)
        {
            const auto size = x.ids.size();
            const auto tmp = x.fn(state, p, w);
            errorif(tmp.size() != size, "You have provided a system of equation constraints whose size, ", tmp.size(), ", does not match with the number of constraint ids, ", size, ", whose ids being: ", x.ids);
            res.segment(offset, size) = tmp;
            offset += size;
        }

        return res;
    };

    return econstraints;
}

/// Assemble the complete system of reactivity constraints to be satisfied at chemical equilibrium.
auto assembleReactivityConstraints(EquilibriumSpecs const& specs) -> ReactivityConstraints
{
    // Create copies of all provided reactivity constraints to avoid full capture of `specs` in lambda function below
    const auto rconstraints_single = specs.reactivityConstraintsSingle();
    const auto rconstraints_system = specs.reactivityConstraintsSystem();

    // The complete system of reactivity constraints that will be created below.
    ReactivityConstraints rconstraints;

    // Collect the id constraints across the single reactivity constraints
    for(auto const& x : rconstraints_single)
        rconstraints.ids.push_back(x.id);

    // Collect the id constraints across the system of reactivity constraints
    for(auto const& x : rconstraints_system)
        rconstraints.ids.insert(rconstraints.ids.end(), x.ids.begin(), x.ids.end());

    // Auxiliary variables
    const auto Nr = rconstraints.ids.size();         // the total number of reactivity constraints
    const auto Nn = specs.system().species().size(); // the number of species in the chemical system
    const auto Np = specs.numControlVariablesP();    // the number of control variables p

    // Create the final constraint coefficient matrices Kn and Kp
    rconstraints.Kn = zeros(Nr, Nn);
    rconstraints.Kp = zeros(Nr, Np);

    // Define auxiliary offset variable to keep track of the current row in the assembly of Kn and Kp below
    auto offset = 0;

    // Go over the single reactivity constraints...
    for(auto const& x : rconstraints_single)
    {
        errorif(x.Kp.size() != 0 && x.Kp.size() != Np, "The size of vector Kp in the reactivity constraint with id `", x.id, "` is ", x.Kp.size(), " which does not match with the number of p control variables, ", Np, ". Ensure the right number of p control variables have been set in the EquilibriumSpecs object.");
        if(x.Kn.size())
            rconstraints.Kn.row(offset) = x.Kn;
        if(x.Kp.size())
            rconstraints.Kp.row(offset) = x.Kp;
        offset += 1;
    }

    // Go over the system of reactivity constraints...
    for(auto const& x : rconstraints_system)
    {
        errorif(x.Kp.size() != 0 && x.Kp.cols() != Np, "The number of columns in matrix Kp in the system of reactivity constraints with ids `", x.ids, "` is ", x.Kp.cols(), " which does not match with the number of p control variables, ", Np, ". Ensure the right number of p control variables have been set in the EquilibriumSpecs object.");
        const auto size = x.ids.size();
        if(rconstraints.Kn.size())
            rconstraints.Kn.middleRows(offset, size) = x.Kn;
        if(rconstraints.Kp.size())
            rconstraints.Kp.middleRows(offset, size) = x.Kp;
        offset += size;
    }

    return rconstraints;
}

} // namespace

EquilibriumSpecs::EquilibriumSpecs(ChemicalSystem const& system)
: m_system(system)
{
    // By default, start with T and P as unknown *p* control variables.
    addControlVariableP({ "T" });
    addControlVariableP({ "P" });

    // By default, start with all surface areas as known *w* input variables.
    for(auto const& surface : system.surfaces())
        addInput("surfaceArea[" + surface.name() + "]"); // e.g., surfaceArea[AqueousPhase:GaseousPhase], surfaceArea[Calcite]
}

//=================================================================================================
//
// STATIC METHODS TO CREATE PREDEFINED CHEMICAL EQUILIBRIUM SPECIFICATIONS
//
//=================================================================================================

auto EquilibriumSpecs::TP(ChemicalSystem const& system) -> EquilibriumSpecs
{
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    return specs;
}

auto EquilibriumSpecs::HP(ChemicalSystem const& system) -> EquilibriumSpecs
{
    EquilibriumSpecs specs(system);
    specs.enthalpy();
    specs.pressure();
    return specs;
}

auto EquilibriumSpecs::TV(ChemicalSystem const& system) -> EquilibriumSpecs
{
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.volume();
    return specs;
}

auto EquilibriumSpecs::UV(ChemicalSystem const& system) -> EquilibriumSpecs
{
    EquilibriumSpecs specs(system);
    specs.internalEnergy();
    specs.volume();
    return specs;
}

auto EquilibriumSpecs::SP(ChemicalSystem const& system) -> EquilibriumSpecs
{
    EquilibriumSpecs specs(system);
    specs.entropy();
    specs.pressure();
    return specs;
}

auto EquilibriumSpecs::SV(ChemicalSystem const& system) -> EquilibriumSpecs
{
    EquilibriumSpecs specs(system);
    specs.entropy();
    specs.volume();
    return specs;
}

//=================================================================================================
//
// METHODS TO SPECIFY THERMODYNAMIC CONSTRAINTS
//
//=================================================================================================

auto EquilibriumSpecs::temperature() -> void
{
    addInput("T");
    pvars = removefn(pvars, RKT_LAMBDA(x, x.name == "T")); // remove any existing *p* control variable for temperature!
}

auto EquilibriumSpecs::pressure() -> void
{
    addInput("P");
    pvars = removefn(pvars, RKT_LAMBDA(x, x.name == "P")); // remove any existing *p* control variable for pressure!
}

auto EquilibriumSpecs::volume() -> void
{
    EquationConstraint constraint;
    constraint.id = "volume";
    const auto idx = addInput("V");
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().volume() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::internalEnergy() -> void
{
    EquationConstraint constraint;
    constraint.id = "internalEnergy";
    const auto idx = addInput("U");
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().internalEnergy() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::enthalpy() -> void
{
    EquationConstraint constraint;
    constraint.id = "enthalpy";
    const auto idx = addInput("H");
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().enthalpy() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::gibbsEnergy() -> void
{
    EquationConstraint constraint;
    constraint.id = "gibbsEnergy";
    const auto idx = addInput("G");
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().gibbsEnergy() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::helmholtzEnergy() -> void
{
    EquationConstraint constraint;
    constraint.id = "helmholtzEnergy";
    const auto idx = addInput("A");
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().helmholtzEnergy() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::entropy() -> void
{
    EquationConstraint constraint;
    constraint.id = "entropy";
    const auto idx = addInput("S");
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().entropy() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::charge() -> void
{
    EquationConstraint constraint;
    constraint.id = "charge";
    const auto idx = addInput(constraint.id);
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().charge() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::elementAmount(StringOrIndex const& element) -> void
{
    const auto ielement = detail::resolveElementIndex(m_system, element);
    const auto elementsymbol = m_system.element(ielement).symbol();
    EquationConstraint constraint;
    constraint.id = "elementAmount[" + elementsymbol + "]";
    const auto idx = addInput(constraint.id);
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().elementAmount(ielement) - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::elementAmountInPhase(StringOrIndex const& element, StringOrIndex const& phase) -> void
{
    const auto ielement = detail::resolveElementIndex(m_system, element);
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto elementsymbol = m_system.element(ielement).symbol();
    const auto phasename = m_system.phase(iphase).name();
    EquationConstraint constraint;
    constraint.id = "elementAmountInPhase[" + elementsymbol + "][" + phasename + "]";
    const auto idx = addInput(constraint.id);
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().elementAmountInPhase(ielement, iphase) - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::elementMass(StringOrIndex const& element) -> void
{
    const auto ielement = detail::resolveElementIndex(m_system, element);
    const auto elementsymbol = m_system.element(ielement).symbol();
    EquationConstraint constraint;
    constraint.id = "elementMass[" + elementsymbol + "]";
    const auto idx = addInput(constraint.id);
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().elementMass(ielement) - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::elementMassInPhase(StringOrIndex const& element, StringOrIndex const& phase) -> void
{
    const auto ielement = detail::resolveElementIndex(m_system, element);
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto elementsymbol = m_system.element(ielement).symbol();
    const auto phasename = m_system.phase(iphase).name();
    EquationConstraint constraint;
    constraint.id = "elementMassInPhase[" + elementsymbol + "][" + phasename + "]";
    const auto idx = addInput(constraint.id);
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().elementMassInPhase(ielement, iphase) - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::phaseAmount(StringOrIndex const& phase) -> void
{
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto phasename = m_system.phase(iphase).name();
    EquationConstraint constraint;
    constraint.id = "phaseAmount[" + phasename + "]";
    const auto idx = addInput(constraint.id);
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().phaseProps(iphase).amount() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::phaseMass(StringOrIndex const& phase) -> void
{
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto phasename = m_system.phase(iphase).name();
    EquationConstraint constraint;
    constraint.id = "phaseMass[" + phasename + "]";
    const auto idx = addInput(constraint.id);
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().phaseProps(iphase).mass() - w[idx];
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::phaseVolume(StringOrIndex const& phase) -> void
{
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto phasename = m_system.phase(iphase).name();
    EquationConstraint constraint;
    constraint.id = "phaseVolume[" + phasename + "]";
    const auto idx = addInput(constraint.id);
    constraint.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().phaseProps(iphase).volume() - w[idx];
    };
    addConstraint(constraint);
}

//=================================================================================================
//
// METHODS TO SPECIFY SURFACE AREA CONDITIONS
//
//=================================================================================================

auto EquilibriumSpecs::surfaceAreas() -> void
{
    for(auto i = 0; i < system().surfaces().size(); ++i)
        surfaceArea(i);
}

auto EquilibriumSpecs::surfaceArea(StringOrIndex const& surface) -> void
{
    auto const& surfaces = m_system.surfaces();
    auto const isurface = detail::resolveSurfaceIndex(surfaces, surface);
    errorif(isurface == surfaces.size(), "There is no surface area with name or index `", stringfy(surface), "` in the chemical system.");
    auto const id = "surfaceArea[" + surfaces[isurface].name() + "]";
    addInput(id);
    pvars = removefn(pvars, RKT_LAMBDA(x, x.name == id)); // remove any existing *p* control variable for this surface area!
}

//=================================================================================================
//
// METHODS TO SPECIFY UNKNOWN INPUT CONDITIONS
//
//=================================================================================================

auto EquilibriumSpecs::unknownTemperature() -> void
{
    if(!isTemperatureUnknown())
    {
        addControlVariableP({ "T" });
        m_inputs = remove(m_inputs, "T"); // remove the existing input for temperature
    }
}

auto EquilibriumSpecs::unknownPressure() -> void
{
    if(!isPressureUnknown())
    {
        addControlVariableP({ "P" });
        m_inputs = remove(m_inputs, "P"); // remove the existing input for pressure
    }
}

auto EquilibriumSpecs::unknownSurfaceAreas() -> void
{
    for(auto i = 0; i < system().surfaces().size(); ++i)
        unknownSurfaceArea(i);
}

auto EquilibriumSpecs::unknownSurfaceArea(StringOrIndex const& surface) -> void
{
    auto const& surfaces = m_system.surfaces();
    auto const isurface = detail::resolveSurfaceIndex(surfaces, surface);
    errorif(isurface == surfaces.size(), "There is no surface area with name or index `", stringfy(surface), "` in the chemical system.");
    auto const id = "surfaceArea[" + surfaces[isurface].name() + "]";

    if(!isSurfaceAreaUnknown(isurface))
    {
        addControlVariableP({ id });
        m_inputs = remove(m_inputs, id); // remove the existing input for surface area
    }
}

//=================================================================================================
//
// METHODS TO SPECIFY CHEMICAL POTENTIAL CONSTRAINTS
//
//=================================================================================================

auto EquilibriumSpecs::chemicalPotential(String substance) -> void
{
    const auto pid = "u[" + substance + "]";
    const auto idx = addInput(pid);
    ControlVariableQ qvar;
    qvar.name = "[" + substance + "]";
    qvar.substance = substance;
    qvar.id = pid;
    qvar.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return w[idx];
    };
    addControlVariableQ(qvar);
}

auto EquilibriumSpecs::lnActivity(Species const& species) -> void
{
    const auto pid = "ln(a[" + species.name() + "])";
    const auto idx = addInput(pid);
    ControlVariableQ qvar;
    qvar.name = "[" + species.name() + "]";
    qvar.substance = species.formula();
    qvar.id = pid;
    qvar.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        auto const& T = state.temperature();
        auto const& P = state.pressure();
        auto const& R = universalGasConstant;
        const auto u0 = species.standardThermoProps(T, P).G0;
        return u0 + R*T*w[idx];
    };
    addControlVariableQ(qvar);
}

auto EquilibriumSpecs::lnActivity(String name) -> void
{
    const auto specieslist = m_system.database().species();
    const auto idx = specieslist.findWithName(name);
    errorif(idx >= specieslist.size(),
        "Could not impose an activity constraint for species with name `", name, "` "
        "because it is not in the database.");
    const auto species = specieslist[idx];
    lnActivity(species);
}

auto EquilibriumSpecs::lgActivity(String name) -> void
{
    lnActivity(name);
}

auto EquilibriumSpecs::activity(String name) -> void
{
    lnActivity(name);
}

auto EquilibriumSpecs::fugacity(String gas) -> void
{
    const auto species = getGaseousSpecies(m_system.database(), gas);
    errorif(species.name().empty(),
        "Could not impose the fugacity constraint for gas `", gas, "` because "
        "there is no gaseous species in the database with this chemical formula.");
    const auto pid = "f[" + gas + "]";
    const auto idx = addInput(pid);
    ControlVariableQ qvar;
    qvar.name = "[" + species.name() + "]";
    qvar.substance = species.formula();
    qvar.id = pid;
    qvar.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        auto const& T = state.temperature();
        auto const& P = state.pressure();
        auto const& R = universalGasConstant;
        const auto u0 = species.standardThermoProps(T, P).G0;
        return u0 + R*T*log(w[idx]);
    };
    addControlVariableQ(qvar);
}

auto EquilibriumSpecs::pH() -> void
{
    const auto species = getAqueousSpecies(m_system.database(), "H+");
    errorif(species.name().empty(),
        "Could not impose pH constraint because the database has "
        "no aqueous species with chemical formula `H+`.");
    const auto pid = "pH";
    const auto idx = addInput(pid);
    ControlVariableQ qvar;
    qvar.name = "[H+]";
    qvar.substance = "H+";
    qvar.id = pid;
    qvar.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        auto const& T = state.temperature();
        auto const& P = state.pressure();
        auto const& R = universalGasConstant;
        const auto u0 = species.standardThermoProps(T, P).G0;
        const auto pH = w[idx];
        return u0 + R*T*(-pH*ln10);
    };
    addControlVariableQ(qvar);
}

auto EquilibriumSpecs::pMg() -> void
{
    const auto species = getAqueousSpecies(m_system.database(), "Mg+2");
    errorif(species.name().empty(),
        "Could not impose pMg constraint because the database has "
        "no aqueous species with chemical formula `Mg+2`.");
    const auto pid = "pMg";
    const auto idx = addInput(pid);
    ControlVariableQ qvar;
    qvar.name = "[Mg+2]";
    qvar.substance = "Mg+2";
    qvar.id = pid;
    qvar.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        auto const& T  = state.temperature();
        auto const& P  = state.pressure();
        auto const& R  = universalGasConstant;
        const auto u0  = species.standardThermoProps(T, P).G0;
        const auto pMg = w[idx];
        return u0 + R*T*(-pMg*ln10);
    };
    addControlVariableQ(qvar);
}

auto EquilibriumSpecs::pE() -> void
{
    const auto idx = addInput("pE");
    ControlVariableQ qvar;
    qvar.name = "[e-]";
    qvar.substance = ChemicalFormula("e-");
    qvar.id = "pE";
    qvar.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        auto const& T = state.temperature();
        auto const& R = universalGasConstant;
        const auto pE = w[idx];
        return R*T*(-pE*ln10);
    };
    addControlVariableQ(qvar);
}

auto EquilibriumSpecs::Eh() -> void
{
    const auto F = faradayConstant; // in C/mol
    const auto idx = addInput("Eh");
    ControlVariableQ qvar;
    qvar.name = "[e-]";
    qvar.substance = ChemicalFormula("e-");
    qvar.id = "Eh";
    qvar.fn = [=](ChemicalState const& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return -F * w[idx]; // in J/mol (chemical potential of electron)
    };
    addControlVariableQ(qvar);
}

//=================================================================================================
//
// METHODS TO SPECIFY HOW THE CHEMICAL SYSTEM IS OPEN
//
//=================================================================================================

auto EquilibriumSpecs::openTo(ChemicalFormula const& substance) -> void
{
    addUnknownTitrantAmount(substance);
}

//=================================================================================================
//
// METHODS TO SPECIFY ADDITIONAL UNKNOWNS
//
//=================================================================================================

auto EquilibriumSpecs::addUnknownTitrantAmount(ChemicalFormula const& substance) -> void
{
    ControlVariableP pvar;
    pvar.name = "[" + substance.str() + "]";
    pvar.substance = substance;
    addControlVariableP(pvar);
}

auto EquilibriumSpecs::addUnknownChemicalPotential(String const& species) -> void
{
    errorif(contains(species_with_unknown_chemical_potentials, species), "Cannot add an unknown for the chemical potential of ", species, " because this has already been done directly or indirectly.");
    species_with_unknown_chemical_potentials.push_back(species);
    const auto ispecies = m_system.species().index(species);
    ControlVariableP pvar;
    pvar.name = "u[" + species + "]";
    pvar.ispecies = ispecies;
    pvar.fn = [=](ChemicalState const& state, real const& pk) -> real
    {
        return pk;
    };
    addControlVariableP(pvar);
}

auto EquilibriumSpecs::addUnknownStandardChemicalPotential(String const& species) -> void
{
    errorif(contains(species_with_unknown_chemical_potentials, species), "Cannot add an unknown for the standard chemical potential of ", species, " because this has already been done directly or indirectly.");
    species_with_unknown_chemical_potentials.push_back(species);
    const auto ispecies = m_system.species().index(species);
    ControlVariableP pvar;
    pvar.name = "u0[" + species + "]";
    pvar.ispecies = ispecies;
    pvar.fn = [=](ChemicalState const& state, real const& pk) -> real
    {
        auto const& props = state.props();
        auto const& T = props.temperature();
        auto const& lnai = props.speciesActivityLn(ispecies);
        auto const& R = universalGasConstant;
        return pk + R*T*lnai;
    };
    addControlVariableP(pvar);
}

auto EquilibriumSpecs::addUnknownActivity(String const& species) -> void
{
    errorif(contains(species_with_unknown_chemical_potentials, species), "Cannot add an unknown for the activity of ", species, " because this has already been done directly or indirectly.");
    species_with_unknown_chemical_potentials.push_back(species);
    const auto ispecies = m_system.species().index(species);
    ControlVariableP pvar;
    pvar.name = "a[" + species + "]";
    pvar.ispecies = ispecies;
    pvar.fn = [=](ChemicalState const& state, real const& pk) -> real
    {
        auto const& props = state.props();
        auto const& T = props.temperature();
        auto const& G0 = props.speciesStandardGibbsEnergy(ispecies);
        auto const& R = universalGasConstant;
        return G0 + R*T*log(pk);
    };
    addControlVariableP(pvar);
}

auto EquilibriumSpecs::addUnknownActivityCoefficient(String const& species) -> void
{
    errorif(contains(species_with_unknown_chemical_potentials, species), "Cannot add an unknown for the activity coefficient of ", species, " because this has already been done directly or indirectly.");
    species_with_unknown_chemical_potentials.push_back(species);
    const auto ispecies = m_system.species().index(species);
    ControlVariableP pvar;
    pvar.name = "g[" + species + "]";
    pvar.ispecies = ispecies;
    pvar.fn = [=](ChemicalState const& state, real const& pk) -> real
    {
        auto const& props = state.props();
        auto const& T = props.temperature();
        auto const& G0 = props.speciesStandardGibbsEnergy(ispecies);
        auto const& lnci = props.speciesConcentrationLn(ispecies);
        auto const& R = universalGasConstant;
        return G0 + R*T*(lnci + log(pk));
    };
    addControlVariableP(pvar);
}

//=================================================================================================
//
// METHODS TO GET THE NUMBER OF INTRODUCED CONSTRAINTS, PARAMETERS, AND CONTROL VARIABLES
//
//=================================================================================================

auto EquilibriumSpecs::numInputs() const -> Index
{
    return m_inputs.size();
}

auto EquilibriumSpecs::numParams() const -> Index
{
    return m_params.size();
}

auto EquilibriumSpecs::numControlVariables() const -> Index
{
    return numControlVariablesQ() + numControlVariablesP();
}

auto EquilibriumSpecs::numControlVariablesP() const -> Index
{
    return pvars.size();
}

auto EquilibriumSpecs::numControlVariablesQ() const -> Index
{
    return qvars.size();
}

auto EquilibriumSpecs::numTitrants() const -> Index
{
    return titrants_explicit.size() + titrants_implicit.size();
}

auto EquilibriumSpecs::numTitrantsExplicit() const -> Index
{
    return titrants_explicit.size();
}

auto EquilibriumSpecs::numTitrantsImplicit() const -> Index
{
    return titrants_implicit.size();
}

auto EquilibriumSpecs::numEquationConstraints() const -> Index
{
    return econstraints_ids.size();
}

auto EquilibriumSpecs::numReactivityConstraints() const -> Index
{
    return rconstraints_ids.size();
}

auto EquilibriumSpecs::numConstraints() const -> Index
{
    return numEquationConstraints() + numReactivityConstraints() + numControlVariablesQ();
}

//=================================================================================================
//
// METHODS TO GET THE NAMES OF INTRODUCED CONSTRAINTS, PARAMETERS, AND CONTROL VARIABLES
//
//=================================================================================================

auto EquilibriumSpecs::namesInputs() const -> Strings
{
    return m_inputs;
}

auto EquilibriumSpecs::namesParams() const -> Strings
{
    return vectorize(m_params, RKT_LAMBDA(x, x.id()));
}

auto EquilibriumSpecs::namesControlVariables() const -> Strings
{
    return concatenate(namesControlVariablesP(), namesControlVariablesQ());
}

auto EquilibriumSpecs::namesControlVariablesP() const -> Strings
{
    return vectorize(pvars, RKT_LAMBDA(x, x.name));
}

auto EquilibriumSpecs::namesControlVariablesQ() const -> Strings
{
    return vectorize(qvars, RKT_LAMBDA(x, x.name));
}

auto EquilibriumSpecs::namesTitrants() const -> Strings
{
    return concatenate(namesTitrantsExplicit(), namesTitrantsImplicit());
}

auto EquilibriumSpecs::namesTitrantsExplicit() const -> Strings
{
    return vectorize(titrants_explicit, RKT_LAMBDA(x, "[" + x.str() + "]"));
}

auto EquilibriumSpecs::namesTitrantsImplicit() const -> Strings
{
    return vectorize(titrants_implicit, RKT_LAMBDA(x, "[" + x.str() + "]"));
}

auto EquilibriumSpecs::namesConstraints() const -> Strings
{
    Strings names = econstraints_ids;
    names = concatenate(names, vectorize(qvars, RKT_LAMBDA(x, x.id)));
    names = concatenate(names, rconstraints_ids);
    return names;
}

//=================================================================================================
//
// METHODS TO ADD CONSTRAINTS AND PARAMETERS
//
//=================================================================================================

auto EquilibriumSpecs::addControlVariableQ(ControlVariableQ const& qvar) -> void
{
    errorif(qvar.name.empty(), "Could not add *q* control variable because its name is empty.");
    errorif(containsfn(qvars, RKT_LAMBDA(x, x.name == qvar.name)), "Could not add *q* control variable with name `", qvar.name, "` because another *q* control variable has already been added with same name.");
    errorif(qvar.fn == nullptr, "Could not add *q* control variable with name `", qvar.name, "` because its chemical potential function is empty.");
    throwErrorIfTitrantHasBeenRegistered(qvar.substance);
    titrants_implicit.push_back(qvar.substance);
    qvars.push_back(qvar);
}

auto EquilibriumSpecs::addControlVariableP(ControlVariableP const& pvar) -> void
{
    const auto Nn = m_system.species().size();
    errorif(pvar.name.empty(), "Could not add *p* control variable because its name is empty.");
    errorif(containsfn(pvars, RKT_LAMBDA(x, x.name == pvar.name)), "Could not add *p* control variable with name `", pvar.name, "` because another *p* control variable has already been added with same name.");
    errorif(pvar.ispecies < Nn && pvar.fn == nullptr, "Could not add *p* control variable with name `", pvar.name, "` because the chemical potential function is missing.");
    errorif(pvar.ispecies >= Nn && pvar.fn, "Could not add *p* control variable with name `", pvar.name, "` because the index of the species whose chemical potential is unknown is missing (index value is ", pvar.ispecies, ", but number of species is ", Nn, ").");
    if(pvar.substance.str().size())
    {
        throwErrorIfTitrantHasBeenRegistered(pvar.substance);
        titrants_explicit.push_back(pvar.substance);
    }
    pvars.push_back(pvar);
}

auto EquilibriumSpecs::addConstraint(EquationConstraint const& constraint) -> void
{
    errorif(contains(econstraints_ids, constraint.id), "Cannot impose a new equation constraint with repeating id (", constraint.id, ").");
    errorif(constraint.id.empty(), "An equation constraint cannot have an empty id.");
    errorif(!constraint.fn, "The equation constraint with id `", constraint.id, "` should not have an empty function.");
    econstraints_ids.push_back(constraint.id);
    econstraints_single.push_back(constraint);
}

auto EquilibriumSpecs::addConstraints(EquationConstraints const& constraints) -> void
{
    for(auto const& constraintid : constraints.ids)
    {
        errorif(contains(econstraints_ids, constraintid), "Cannot impose a new equation constraint with repeating id (", constraintid, ").");
        errorif(constraintid.empty(), "An equation constraint cannot have an empty id.");
        econstraints_ids.push_back(constraintid);
    }
    errorif(!constraints.fn, "The system of equation constraints with ids `", constraints.ids, "` should not have an empty function.");
    econstraints_system.push_back(constraints);
}

auto EquilibriumSpecs::addReactivityConstraint(ReactivityConstraint const& constraint) -> void
{
    const auto Nn = m_system.species().size();
    const auto Np = numControlVariablesP();
    errorif(contains(rconstraints_ids, constraint.id), "Cannot impose a new reactivity constraint with repeating id (", constraint.id, ").");
    errorif(constraint.id.empty(), "A reactivity constraint cannot have an empty id.");
    errorif(constraint.Kn.size() == 0 && constraint.Kp.size() == 0, "The vectors Kn and Kp in the reactivity constraint with id `", constraint.id, "` should not be both empty.");
    errorif(constraint.Kn.size() != 0 && constraint.Kn.size() != Nn, "The size of vector Kn in the reactivity constraint with id `", constraint.id, "` is ", constraint.Kn.size(), " which does not match with the number of species, ", Nn, ".");
    errorif(constraint.Kp.size() != 0 && constraint.Kp.size() != Nn, "The size of vector Kp in the reactivity constraint with id `", constraint.id, "` is ", constraint.Kp.size(), " which does not match with the number of p control variables, ", Np, ". Ensure the right number of p control variables have been set first before specifying reactivity constraints in the EquilibriumSpecs object.");
    rconstraints_ids.push_back(constraint.id);
    rconstraints_single.push_back(constraint);
}

auto EquilibriumSpecs::addReactivityConstraints(ReactivityConstraints const& constraints) -> void
{
    const auto Nn = m_system.species().size();
    const auto Np = numControlVariablesP();
    errorif(constraints.Kn.size() == 0 && constraints.Kp.size() == 0, "The matrices Kn and Kp in the system of reactivity constraints with ids `", constraints.ids, "` should not be both empty.");
    errorif(constraints.Kn.size() != 0 && constraints.Kn.cols() != Nn, "The number of columns in matrix Kn in the system of reactivity constraints with ids `", constraints.ids, "` is ", constraints.Kn.cols(), " which does not match with the number of species, ", Nn, ".");
    errorif(constraints.Kp.size() != 0 && constraints.Kp.cols() != Np, "The number of columns in matrix Kp in the system of reactivity constraints with ids `", constraints.ids, "` is ", constraints.Kp.cols(), " which does not match with the number of p control variables, ", Np, ". Ensure the right number of p control variables have been set first before specifying reactivity constraints in the EquilibriumSpecs object.");
    for(auto const& constraintid : constraints.ids)
    {
        errorif(contains(rconstraints_ids, constraintid), "Cannot impose a new reactivity constraint with repeating id (", constraintid, ").");
        errorif(constraintid.empty(), "A reactivity constraint cannot have an empty id.");
        rconstraints_ids.push_back(constraintid);
    }
    errorif(constraints.Kn.size() != 0 && constraints.Kn.rows() != constraints.ids.size(), "The number of rows in matrix Kn in the system of reactivity constraints with ids `", constraints.ids, "` is ", constraints.Kn.rows(), " which does not match with the number of constraint ids, ", constraints.ids.size(), ".");
    errorif(constraints.Kp.size() != 0 && constraints.Kp.rows() != constraints.ids.size(), "The number of rows in matrix Kp in the system of reactivity constraints with ids `", constraints.ids, "` is ", constraints.Kp.rows(), " which does not match with the number of constraint ids, ", constraints.ids.size(), ".");
    rconstraints_system.push_back(constraints);
}

auto EquilibriumSpecs::addInput(String const& var) -> Index
{
    errorif(contains(m_inputs, var), "Could not add input variable with id `", var, "` to the "
        "list of input variables for the equilibrium problem because "
        "another input already exists with same id.");
    m_inputs.push_back(var);
    return m_inputs.size() - 1; // the index of the just added input variable
}

auto EquilibriumSpecs::addInput(Param const& param) -> Index
{
    const auto idx = addInput(param.id());
    m_params.push_back(param);
    m_params_idxs.push_back(idx);
    return idx;
}

//=================================================================================================
//
// MISCELLANEOUS METHODS
//
//=================================================================================================

auto EquilibriumSpecs::system() const -> ChemicalSystem const&
{
    return m_system;
}

auto EquilibriumSpecs::inputs() const -> Strings const&
{
    return m_inputs;
}

auto EquilibriumSpecs::params() const -> const Vec<Param>&
{
    return m_params;
}

auto EquilibriumSpecs::indicesParams() const -> const Vec<Index>&
{
    return m_params_idxs;
}

auto EquilibriumSpecs::isTemperatureUnknown() const -> bool
{
    return containsfn(pvars, RKT_LAMBDA(x, x.name == "T"));
}

auto EquilibriumSpecs::isPressureUnknown() const -> bool
{
    return containsfn(pvars, RKT_LAMBDA(x, x.name == "P"));
}

auto EquilibriumSpecs::isSurfaceAreaUnknown(StringOrIndex const& surface) const -> bool
{
    auto const& surfaces = m_system.surfaces();
    auto const isurface = detail::resolveSurfaceIndex(surfaces, surface);
    errorif(isurface == surfaces.size(), "There is no surface area with name or index `", stringfy(surface), "` in the chemical system.");
    auto const id = "surfaceArea[" + surfaces[isurface].name() + "]";
    return containsfn(pvars, RKT_LAMBDA(x, x.name == id));
}

auto EquilibriumSpecs::indexTemperatureAmongInputVariables() const -> Index
{
    auto k = index(m_inputs, "T");
    return k < m_inputs.size() ? k : Index(-1);
}

auto EquilibriumSpecs::indexTemperatureAmongControlVariablesP() const -> Index
{
    auto k = indexfn(pvars, RKT_LAMBDA(x, x.name == "T"));
    return k < pvars.size() ? k : Index(-1);
}

auto EquilibriumSpecs::indexPressureAmongInputVariables() const -> Index
{
    auto k = index(m_inputs, "P");
    return k < m_inputs.size() ? k : Index(-1);
}

auto EquilibriumSpecs::indexPressureAmongControlVariablesP() const -> Index
{
    auto k = indexfn(pvars, RKT_LAMBDA(x, x.name == "P"));
    return k < pvars.size() ? k : Index(-1);
}

auto EquilibriumSpecs::indicesSurfaceAreasAmongInputVariables() const -> Indices
{
    Indices indices;
    for(auto const& surface : system().surfaces())
    {
        auto id = "surfaceArea[" + surface.name() + "]";
        auto idx = index(m_inputs, "surfaceArea[" + surface.name() + "]");
        if(idx < m_inputs.size())
            indices.push_back(idx);
    }
    return indices;
}

auto EquilibriumSpecs::indicesSurfaceAreasAmongControlVariablesP() const -> Indices
{
    Indices indices;
    for(auto const& surface : system().surfaces())
    {
        auto id = "surfaceArea[" + surface.name() + "]";
        auto idx = indexfn(pvars, RKT_LAMBDA(x, x.name == id));
        if(idx < m_inputs.size())
            indices.push_back(idx);
    }
    return indices;
}

auto EquilibriumSpecs::indicesSurfaceAreasKnown() const -> Indices
{
    Indices indices;
    for(auto const& [i, surface] : enumerate(system().surfaces()))
        if(contains(m_inputs, "surfaceArea[" + surface.name() + "]"))
            indices.push_back(i);
    return indices;
}

auto EquilibriumSpecs::indicesSurfaceAreasUnknown() const -> Indices
{
    Indices indices;
    for(auto const& [i, surface] : enumerate(system().surfaces()))
    {
        auto id = "surfaceArea[" + surface.name() + "]";
        if(containsfn(pvars, RKT_LAMBDA(x, x.name == id)))
            indices.push_back(i);
    }
    return indices;
}

auto EquilibriumSpecs::indexInputVariable(String const& name) const -> Index
{
    return index(m_inputs, name);
}

auto EquilibriumSpecs::indexControlVariableP(String const& name) const -> Index
{
    return indexfn(pvars, RKT_LAMBDA(x, x.name == name));
}

auto EquilibriumSpecs::indexControlVariableQ(String const& name) const -> Index
{
    return indexfn(qvars, RKT_LAMBDA(x, x.name == name));
}

auto EquilibriumSpecs::controlVariablesQ() const -> const Vec<ControlVariableQ>&
{
    return qvars;
}

auto EquilibriumSpecs::controlVariablesP() const -> const Vec<ControlVariableP>&
{
    return pvars;
}

auto EquilibriumSpecs::titrants() const -> Vec<ChemicalFormula>
{
    return concatenate(titrants_explicit, titrants_implicit);
}

auto EquilibriumSpecs::titrantsExplicit() const -> Vec<ChemicalFormula>
{
    return titrants_explicit;
}

auto EquilibriumSpecs::titrantsImplicit() const -> Vec<ChemicalFormula>
{
    return titrants_implicit;
}

auto EquilibriumSpecs::equationConstraintsSingle() const -> Vec<EquationConstraint> const&
{
    return econstraints_single;
}

auto EquilibriumSpecs::equationConstraintsSystem() const -> Vec<EquationConstraints> const&
{
    return econstraints_system;
}

auto EquilibriumSpecs::equationConstraints() const -> EquationConstraints
{
    return assembleEquationConstraints(*this);
}

auto EquilibriumSpecs::reactivityConstraintsSingle() const -> Vec<ReactivityConstraint> const&
{
    return rconstraints_single;
}

auto EquilibriumSpecs::reactivityConstraintsSystem() const -> Vec<ReactivityConstraints> const&
{
    return rconstraints_system;
}

auto EquilibriumSpecs::reactivityConstraints() const -> ReactivityConstraints
{
    return assembleReactivityConstraints(*this);
}

//=================================================================================================
//
// PRIVATE METHODS
//
//=================================================================================================

auto EquilibriumSpecs::throwErrorIfTitrantHasBeenRegistered(ChemicalFormula const& substance) const -> void
{
    const auto already_an_explicit_titrant = containsfn(titrants_explicit, RKT_LAMBDA(x, x.equivalent(substance)));
    errorif(already_an_explicit_titrant, "Cannot specify again that the chemical system is open to substance ", substance.str(), ". "
        "You have already explicitly specified that the chemical system is open to this substance.");
    const auto idx = indexfn(titrants_implicit, RKT_LAMBDA(x, x.equivalent(substance)));
    const auto already_an_implicit_titrant = idx < titrants_implicit.size();
    errorif(already_an_implicit_titrant, "Cannot specify again that the chemical system is open to substance ", substance.str(), ". "
        "You have implicitly specified that the chemical system is open to this substance when imposing a chemical potential constraint for it.");
}

} // namespace Reaktoro
