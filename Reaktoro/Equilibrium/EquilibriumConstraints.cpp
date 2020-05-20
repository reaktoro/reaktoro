// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "EquilibriumConstraints.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionEquation.hpp>

namespace Reaktoro {

//=================================================================================================
//
// EquilibriumConstraints::Impl
//
//=================================================================================================

struct EquilibriumConstraints::Impl
{
    // ConstraintFn gT;
    // ConstraintFn gP;

    // /// The derivatives of the objective function f with respect to x = n.
    // Vec<ConstraintFn> dfdx;

    // /// The derivatives of the objective function f with respect to p = (T, P).
    // Vec<ConstraintFn> dfdp;

    // /// The derivatives of the objective function f with respect to q (the amounts of titrants).
    // Vec<ConstraintFn> dfdq;

    /// The chemical system associated with the equilibrium constraints.
    ChemicalSystem system;

    /// The list of general equilibrium constraints to be imposed.
    Vec<EquilibriumConstraint> eqconstraints;

    /// The list of chemical potential constraints to be imposed.
    Vec<ChemicalPotentialConstraint> uconstraints;

    /// The reactivity constraints to be imposed.
    ReactivityConstraints rconstraints;

    /// Construct a EquilibriumConstraints::Impl object.
    Impl(const ChemicalSystem& system)
    : system(system)
    {}

    /// Return a Control object to initiate the imposition of a general equilibrium constraint.
    auto control() -> Control
    {
        return Control(eqconstraints);
    }

    /// Return a Fix object to initiate the imposition of a chemical potential constraint.
    auto fix() -> Fix
    {
        uconstraints.push_back(ChemicalPotentialConstraint());
        return Fix(system, uconstraints.back());
    }

    /// Return a Prevent object to initiate the imposition of a reactivity constraint.
    auto prevent() -> Prevent
    {
        return Prevent(system, rconstraints);
    }

    /// Return the conservation matrix associated with the equilibrium constraints.
    auto conservationMatrix() -> MatrixXd
    {
        // const auto num_elements = system.elements().size();
        // const auto num_species = system.species().size();
        // const auto num_inert_reactions = rconstraints.reactions_cannot_react.size();
        // const auto num_charge = 1;
        // const auto num_titrants = eqconstraints.size() + uconstraints.size();

        // const auto num_rows = num_elements + num_charge + num_inert_reactions;
        // const auto num_columns = num_species + num_titrants;

        // const auto& A = system.formulaMatrix();
        // const auto nrows = A.rows() + rconstraints.reactions_cannot_react.size();
        // MatrixXd C(nrows, ncols);
    }
};

//=================================================================================================
//
// EquilibriumConstraints
//
//=================================================================================================

EquilibriumConstraints::EquilibriumConstraints(const ChemicalSystem& system)
{

}

auto EquilibriumConstraints::control() -> Control
{

}

auto EquilibriumConstraints::fix() -> Fix
{

}

auto EquilibriumConstraints::prevent() -> Prevent
{

}

//=================================================================================================
//
// EquilibriumConstraints::Control
//
//=================================================================================================

EquilibriumConstraints::Control::Control(Vec<EquilibriumConstraint>& constraints)
: constraints(constraints)
{}

auto EquilibriumConstraints::Control::temperature() -> Until
{
    constraints.push_back(EquilibriumConstraint());
    constraints.back().control = "Temperature";
    return EquilibriumConstraints::Until(constraints.back().fn);
}

auto EquilibriumConstraints::Control::pressure() -> Until
{
    constraints.push_back(EquilibriumConstraint());
    constraints.back().control = "Pressure";
    return EquilibriumConstraints::Until(constraints.back().fn);
}

auto EquilibriumConstraints::Control::titrationOf(String titrant) -> Until
{
    constraints.push_back(EquilibriumConstraint());
    constraints.back().control = titrant;
    return EquilibriumConstraints::Until(constraints.back().fn);
}

auto EquilibriumConstraints::Control::titrationOfEither(String titrant1, String titrant2) -> Until
{
    error(true, "Method EquilibriumConstraints::Control::titrationOfEither has not been implemented yet.");
    return EquilibriumConstraints::Until(constraints.back().fn);
}

//=================================================================================================
//
// EquilibriumConstraints::Until
//
//=================================================================================================

EquilibriumConstraints::Until::Until(ConstraintFn& constraintfn)
: constraintfn(constraintfn)
{}

auto EquilibriumConstraints::Until::until() const -> EquilibriumConstraints::Attained
{
    return EquilibriumConstraints::Attained(constraintfn);
}

auto EquilibriumConstraints::Until::until(const ConstraintFn& fn) -> void
{
    constraintfn = fn;
}

//=================================================================================================
//
// EquilibriumConstraints::Attained
//
//=================================================================================================

EquilibriumConstraints::Attained::Attained(ConstraintFn& constraintfn)
: constraintfn(constraintfn)
{}

auto EquilibriumConstraints::Attained::volume(real value, String unit) -> void
{
    value = units::convert(value, unit, "m3");
    constraintfn = [=](const ChemicalProps& props) { return props.volume() - value; };
}

auto EquilibriumConstraints::Attained::internalEnergy(real value, String unit) -> void
{
    value = units::convert(value, unit, "J");
    constraintfn = [=](const ChemicalProps& props) { return props.internalEnergy() - value; };
}

auto EquilibriumConstraints::Attained::enthalpy(real value, String unit) -> void
{
    value = units::convert(value, unit, "J");
    constraintfn = [=](const ChemicalProps& props) { return props.enthalpy() - value; };
}

auto EquilibriumConstraints::Attained::gibbsEnergy(real value, String unit) -> void
{
    value = units::convert(value, unit, "J");
    constraintfn = [=](const ChemicalProps& props) { return props.gibbsEnergy() - value; };
}

auto EquilibriumConstraints::Attained::helmholtzEnergy(real value, String unit) -> void
{
    value = units::convert(value, unit, "J");
    constraintfn = [=](const ChemicalProps& props) { return props.helmholtzEnergy() - value; };
}

auto EquilibriumConstraints::Attained::entropy(real value, String unit) -> void
{
    value = units::convert(value, unit, "J/K");
    constraintfn = [=](const ChemicalProps& props) { return props.entropy() - value; };
}

//=================================================================================================
//
// EquilibriumConstraints::Fix
//
//=================================================================================================

namespace detail {

/// Create a chemical potential constraint with given species and its constrained ln activity.
/// @param species
auto lnActivityConstraint(const Species& species, real value) -> ChemicalPotentialConstraint
{
    const auto R = universalGasConstant;
    ChemicalPotentialConstraint constraint;
    constraint.formula = species.formula();
    constraint.fn = [=](const ChemicalProps& props) {
        const auto T = props.temperature();
        const auto P = props.pressure();
        const auto u0 = species.props(T, P).G0; // the standard chemical potential of the species
        return u0 + R*T*value;
    };
    return constraint;
}

} // namespace detail

EquilibriumConstraints::Fix::Fix(const ChemicalSystem& system, ChemicalPotentialConstraint& uconstraint)
: system(system), uconstraint(uconstraint)
{}

auto EquilibriumConstraints::Fix::chemicalPotential(String substance, real value, String unit) -> void
{
    value = units::convert(value, unit, "J/mol");
    uconstraint.formula = ChemicalFormula(substance);
    uconstraint.fn = [=](const ChemicalProps& props) { return value; };
}

auto EquilibriumConstraints::Fix::lnActivity(String name, real value) -> void
{
    const auto idx = system.database().species().findWithName(name);
    error(idx >= system.database().species().size(),
        "Cannot impose an activity constraint for species `", name, "` "
        "because it is not in the database.");
    const auto species = system.database().species()[idx];
    uconstraint = detail::lnActivityConstraint(species, value);
}

auto EquilibriumConstraints::Fix::lgActivity(String name, real value) -> void
{
    lnActivity(name, value * ln10);
}

auto EquilibriumConstraints::Fix::activity(String name, real value) -> void
{
    lnActivity(name, log(value));
}

auto EquilibriumConstraints::Fix::fugacity(String name, real value, String unit) -> void
{
    value = units::convert(value, unit, "bar");
    const auto& gases = system.database().speciesWithAggregateState(AggregateState::Gas);
    const auto idx = gases.findWithName(name);
    error(idx >= gases.size(),
        "Cannot impose a fugacity constraint for species `", name, "` "
        "because it is not a gaseous species in the database.");
    const auto species = gases[idx];
    uconstraint = detail::lnActivityConstraint(species, value);
}

auto EquilibriumConstraints::Fix::pH(real value) -> void
{
    const auto aqspecies = system.database().speciesWithAggregateState(AggregateState::Aqueous);
    const auto idx = aqspecies.findWithFormula("H+");
    error(idx >= aqspecies.size(),
        "Cannot impose pH constraint because the database has "
        "no aqueous species with chemical formula `H+`.");
    const auto species = aqspecies[idx];
    uconstraint = detail::lnActivityConstraint(species, -value * ln10); // pH = -log10(a[H+]) => ln(a[H+]) = -pH * ln10
}

auto EquilibriumConstraints::Fix::pMg(real value) -> void
{
    const auto aqspecies = system.database().speciesWithAggregateState(AggregateState::Aqueous);
    const auto idx = aqspecies.findWithFormula("Mg+2");
    error(idx >= aqspecies.size(),
        "Cannot impose pMg constraint because the database has "
        "no aqueous species with chemical formula `Mg+2`.");
    const auto species = aqspecies[idx];
    uconstraint = detail::lnActivityConstraint(species, -value * ln10); // pMg = -log10(a[Mg+2]) => ln(a[Mg+2]) = -pH * ln10
}

auto EquilibriumConstraints::Fix::pe(real value) -> void
{
    const auto ue0 = 0.0;            // the standard chemical potential of the electron species
    const auto lnae = -value * ln10; // the ln activity of the electron species
    const auto R = universalGasConstant;
    uconstraint.formula = ChemicalFormula("e-");
    uconstraint.fn = [=](const ChemicalProps& props) {
        const auto T = props.temperature();
        return ue0 + R*T*lnae;
    };
}

auto EquilibriumConstraints::Fix::Eh(real value, String unit) -> void
{
    value = units::convert(value, unit, "V"); // in V = J/C
    const auto F = faradayConstant; // in C/mol
    const auto ue = -F * value; // in J/mol (chemical potential of electron)
    uconstraint.formula = ChemicalFormula("e-");
    uconstraint.fn = [=](const ChemicalProps& props) { return ue; };
}

//=================================================================================================
//
// EquilibriumConstraints::Prevent
//
//=================================================================================================

EquilibriumConstraints::Prevent::Prevent(const ChemicalSystem& system, ReactivityConstraints& constraints)
: system(system), constraints(constraints)
{}

auto EquilibriumConstraints::Prevent::fromReacting(Index ispecies) -> void
{
    const auto size = system.species().size();
    error(ispecies >= size, "The given species index ", ispecies, " is "
        "out of bounds, since there are only ", size, " species in the chemical system.");
    constraints.species_cannot_react.insert(ispecies);
}

auto EquilibriumConstraints::Prevent::fromReacting(Pairs<Index, double> equation) -> void
{
    constraints.reactions_cannot_react.push_back(equation);
}

auto EquilibriumConstraints::Prevent::fromReacting(String what) -> void
{
    if(what.find("=") == String::npos) {
        fromReacting(system.species().indexWithName(what));
    }
    else {
        Pairs<Index, double> pairs;
        for(auto [name, coeff] : parseReactionEquation(what))
            pairs.emplace_back(system.species().indexWithName(name), coeff);
        fromReacting(pairs);
    }
}

auto EquilibriumConstraints::Prevent::fromIncreasing(Index ispecies) -> void
{
    const auto size = system.species().size();
    error(ispecies >= size, "The given species index ", ispecies, " is "
        "out of bounds, since there are only ", size, " species in the chemical system.");
    constraints.species_cannot_increase.insert(ispecies);
}

auto EquilibriumConstraints::Prevent::fromIncreasing(String species) -> void
{
    fromIncreasing(system.species().indexWithName(species));
}

auto EquilibriumConstraints::Prevent::fromDecreasing(Index ispecies) -> void
{
    const auto size = system.species().size();
    error(ispecies >= size, "The given species index ", ispecies, " is "
        "out of bounds, since there are only ", size, " species in the chemical system.");
    constraints.species_cannot_decrease.insert(ispecies);
}

auto EquilibriumConstraints::Prevent::fromDecreasing(String species) -> void
{
    fromDecreasing(system.species().indexWithName(species));
}

} // namespace Reaktoro
