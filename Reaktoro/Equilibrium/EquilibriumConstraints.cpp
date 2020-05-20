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
    /// The chemical system associated with the equilibrium constraints.
    ChemicalSystem system;

    /// The imposed equilibrium constraints.
    Data data;

    /// Construct a EquilibriumConstraints::Impl object.
    Impl(const ChemicalSystem& system)
    : system(system)
    {}

    /// Return a Control object to initiate the imposition of a general equilibrium constraint.
    auto control() -> Control
    {
        return Control(data.fconstraints);
    }

    /// Return a Fix object to initiate the imposition of a chemical potential constraint.
    auto fix() -> Fix
    {
        data.uconstraints.push_back({});
        return Fix(system, data.uconstraints.back());
    }

    /// Return a Prevent object to initiate the imposition of a reactivity constraint.
    auto prevent() -> Prevent
    {
        return Prevent(system, data.rconstraints);
    }

    /// Assemble the matrix block A in the conservation matrix C.
    auto assembleMatrixA(MatrixXdRef A) -> void
    {
        A = system.formulaMatrix();
    }

    /// Assemble the matrix block B = [BT BP Bq] in the conservation matrix C.
    auto assembleMatrixB(MatrixXdRef B) -> void
    {
        const auto num_elements = system.elements().size();
        assert(B.cols() == 2 + data.fconstraints.dfdq.size()); // T, P and number of titrants
        assert(B.rows() == 1 + num_elements);              // number of elements plus charge

        auto fill_matrix_col = [&](const auto& formula, auto col)
        {
            col[num_elements] = formula.charge(); // last entry in the column vector is charge of substance
            for(const auto& [element, coeff] : formula.elements()) {
                const auto ielem = system.elements().index(element);
                col[ielem] = coeff;
            }
        };

        auto j = 2; // skip columns BT and BP, since these are zeros
        for(const auto& [formula, fn] : data.fconstraints.dfdq)
            fill_matrix_col(formula, B.col(j++));
    }

    /// Assemble the matrix block S in the conservation matrix C.
    auto assembleMatrixS(MatrixXdRef S) -> void
    {
        assert(S.rows() == data.rconstraints.reactions_cannot_react.size());

        auto fill_matrix_row = [&](const auto& pairs, auto row)
        {
            for(auto [ispecies, coeff] : pairs)
                row[ispecies] = coeff;
        };

        auto i = 0;
        for(const auto& pairs : data.rconstraints.reactions_cannot_react)
            fill_matrix_row(pairs, S.row(i++));
    }

    /// Return the conservation matrix associated with the equilibrium constraints.
    /// The conservation matrix is:
    ///
    /// C = [ A B ]
    ///     [ S 0 ]
    ///
    /// where A is the formula matrix of the species with respect to elements
    /// and charge; B = [BT BP Bq], with BT and BP being zero column vectors
    /// and Bq the formula matrix of the introduced titrants whose amounts are
    /// controlled to attain imposed equilibrium constraints; and S is the
    /// stoichiometric matrix of reactions that cannot progress during the
    /// equilibrium calculation (inert reactions).
    auto conservationMatrix() -> MatrixXd
    {
        const auto num_elements = system.elements().size();
        const auto num_species = system.species().size();
        const auto num_charge = 1;
        const auto num_inert_reactions = data.rconstraints.reactions_cannot_react.size();
        const auto num_controls = 2 + data.fconstraints.dfdq.size();

        const auto num_rows = num_elements + num_charge + num_inert_reactions;
        const auto num_cols = num_species + num_controls;

        MatrixXd C = MatrixXd::Zero(num_rows, num_cols);

        auto A = C.topLeftCorner(num_elements + num_charge, num_species);
        auto B = C.topRightCorner(num_elements + num_charge, num_controls);
        auto S = C.bottomLeftCorner(num_inert_reactions, num_species);

        assembleMatrixA(A);
        assembleMatrixB(B);
        assembleMatrixS(S);

        return C;
    }
};

//=================================================================================================
//
// EquilibriumConstraints
//
//=================================================================================================

EquilibriumConstraints::EquilibriumConstraints(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumConstraints::~EquilibriumConstraints()
{}

auto EquilibriumConstraints::control() -> Control
{
    return pimpl->control();
}

auto EquilibriumConstraints::fix() -> Fix
{
    return pimpl->fix();
}

auto EquilibriumConstraints::prevent() -> Prevent
{
    return pimpl->prevent();
}

auto EquilibriumConstraints::data() const -> const Data&
{
    return pimpl->data;
}

//=================================================================================================
//
// EquilibriumConstraints::Control
//
//=================================================================================================

EquilibriumConstraints::Control::Control(FunctionalConstraints& fconstraints)
: fconstraints(fconstraints)
{}

auto EquilibriumConstraints::Control::temperature() -> Until
{
    return EquilibriumConstraints::Until(fconstraints.dfdT);
}

auto EquilibriumConstraints::Control::pressure() -> Until
{
    return EquilibriumConstraints::Until(fconstraints.dfdP);
}

auto EquilibriumConstraints::Control::titrationOf(String titrant) -> Until
{
    // TODO: Adapt Optima lib so that Ax + Bq = b can be imposed by B^T
    // (transpose) is not considered in the first-order conditions for minimum.
    error(true, "Method EquilibriumConstraints::Control::titrationOf is not supported yet.");
    fconstraints.dfdq.push_back({});
    auto& [formula, fn] = fconstraints.dfdq.back();
    formula = ChemicalFormula(titrant);
    return EquilibriumConstraints::Until(fn);
}

auto EquilibriumConstraints::Control::titrationOfEither(String titrant1, String titrant2) -> Until
{
    error(true, "Method EquilibriumConstraints::Control::titrationOfEither has not been implemented yet.");
    EquilibriumConstraintFn fn;
    return EquilibriumConstraints::Until(fn);
}

//=================================================================================================
//
// EquilibriumConstraints::Until
//
//=================================================================================================

EquilibriumConstraints::Until::Until(EquilibriumConstraintFn& constraintfn)
: constraintfn(constraintfn)
{}

auto EquilibriumConstraints::Until::until() const -> EquilibriumConstraints::Attained
{
    return EquilibriumConstraints::Attained(constraintfn);
}

auto EquilibriumConstraints::Until::until(const EquilibriumConstraintFn& customfn) -> void
{
    constraintfn = customfn;
}

//=================================================================================================
//
// EquilibriumConstraints::Attained
//
//=================================================================================================

EquilibriumConstraints::Attained::Attained(EquilibriumConstraintFn& constraintfn)
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
    constraint.fn = [=](real T, real P) {
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
    uconstraint.fn = [=](real T, real P) { return value; };
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
    uconstraint.fn = [=](real T, real P) { return ue0 + R*T*lnae; };
}

auto EquilibriumConstraints::Fix::Eh(real value, String unit) -> void
{
    value = units::convert(value, unit, "V"); // in V = J/C
    const auto F = faradayConstant;           // in C/mol
    const auto ue = -F * value;               // in J/mol (chemical potential of electron)
    uconstraint.formula = ChemicalFormula("e-");
    uconstraint.fn = [=](real T, real P) { return ue; };
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
