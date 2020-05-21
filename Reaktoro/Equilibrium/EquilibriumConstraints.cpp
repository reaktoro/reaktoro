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
        return Control(data);
    }

    /// Return a Until object to initiate the imposition of a functional equilibrium constraint.
    auto until() -> Until
    {
        return Until(data);
    }

    /// Return a Preserve object to initiate the imposition of properties that must be preserved.
    auto preserve() -> Preserve
    {
        return Preserve(data);
    }

    /// Return a Fix object to initiate the imposition of a chemical potential constraint.
    auto fix() -> Fix
    {
        return Fix(system, data);
    }

    /// Return a Prevent object to initiate the imposition of a reactivity constraint.
    auto prevent() -> Prevent
    {
        return Prevent(system, data);
    }

    /// Assemble the vector with the element and charge coefficients of a chemical formula.
    auto assembleFormulaVector(VectorXdRef vec, const ChemicalFormula& formula) -> void
    {
        const auto num_elements = system.elements().size();
        assert(vec.size() == num_elements + 1);
        vec[num_elements] = formula.charge(); // last entry in the column vector is charge of substance
        for(const auto& [element, coeff] : formula.elements()) {
            const auto ielem = system.elements().index(element);
            vec[ielem] = coeff;
        }
    }

    /// Assemble the matrix block A in the conservation matrix C.
    auto assembleMatrixA(MatrixXdRef A) -> void
    {
        A = system.formulaMatrix();
    }

    /// Assemble the matrix block U in the conservation matrix C.
    auto assembleMatrixU(MatrixXdRef U) -> void
    {
        const auto num_elements = system.elements().size();
        const auto num_substances = data.uconstraints.size();

        assert(U.rows() == 1 + num_elements); // number of elements plus charge
        assert(U.cols() == num_substances);   // number of introduced chemical potential constraints

        auto j = 0;
        for(const auto& [formula, _] : data.uconstraints)
            assembleFormulaVector(U.col(j++), formula);
    }

    /// Assemble the matrix block B = [BT BP Bq] in the conservation matrix C.
    auto assembleMatrixB(MatrixXdRef B) -> void
    {
        const auto num_elements = system.elements().size();
        const auto num_controls = data.controls.size();

        assert(B.rows() == 1 + num_elements); // number of elements plus charge
        assert(B.cols() == num_controls);     // number of introduced control variables

        auto j = data.controls.T + data.controls.P; // skip columns BT and BP (if applicable), since these are zeros
        for(const auto& formula : data.controls.titrants)
            assembleFormulaVector(B.col(j++), formula);
    }

    /// Assemble the matrix block S in the conservation matrix C.
    auto assembleMatrixS(MatrixXdRef S) -> void
    {
        auto inert_reactions = data.rconstraints.reactions_cannot_react;

        assert(S.rows() == inert_reactions.size());

        auto fill_matrix_row = [&](const auto& pairs, auto row)
        {
            for(auto [ispecies, coeff] : pairs)
                row[ispecies] = coeff;
        };

        auto i = 0;
        for(const auto& pairs : inert_reactions)
            fill_matrix_row(pairs, S.row(i++));
    }

    /// Return the conservation matrix associated with the equilibrium constraints.
    /// The conservation matrix is:
    ///
    /// C = [ A U B ]
    ///     [ S 0 0 ]
    ///
    /// where A is the formula matrix of the species with respect to elements
    /// and charge; U is the formula matrix of the substances with fixed
    /// chemical potentials; B = [BT BP Bq], with BT and BP being zero column
    /// vectors and Bq the formula matrix of the introduced titrants whose
    /// amounts are controlled to attain imposed equilibrium constraints; and S
    /// is the stoichiometric matrix of reactions that cannot progress during
    /// the equilibrium calculation (inert reactions).
    auto conservationMatrix() -> MatrixXd
    {
        const auto num_elements = system.elements().size();
        const auto num_species = system.species().size();
        const auto num_charge = 1;
        const auto num_inert_reactions = data.rconstraints.reactions_cannot_react.size();
        const auto num_fixed_chemical_potentials = data.uconstraints.size();
        const auto num_controls = data.controls.size();

        const auto num_rows = num_elements + num_charge + num_inert_reactions;
        const auto num_cols = num_species + num_fixed_chemical_potentials + num_controls;

        MatrixXd C = MatrixXd::Zero(num_rows, num_cols);

        auto A = C.topRows(num_elements + num_charge).leftCols(num_species);
        auto U = C.topRows(num_elements + num_charge).middleCols(num_species, num_fixed_chemical_potentials);
        auto B = C.topRows(num_elements + num_charge).rightCols(num_controls);
        auto S = C.bottomLeftCorner(num_inert_reactions, num_species);

        assembleMatrixA(A);
        assembleMatrixU(U);
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

auto EquilibriumConstraints::until() -> Until
{
    return pimpl->until();
}

auto EquilibriumConstraints::preserve() -> Preserve
{
    return pimpl->preserve();
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

EquilibriumConstraints::Control::Control(EquilibriumConstraints::Data& data)
: data(data)
{}

auto EquilibriumConstraints::Control::temperature() -> Control&
{
    data.controls.T = true;
    return *this;
}

auto EquilibriumConstraints::Control::pressure() -> Control&
{
    data.controls.P = true;
    return *this;
}

auto EquilibriumConstraints::Control::titrationOf(String titrant) -> Control&
{
    // TODO: Adapt Optima lib so that Ax + Bq = b can be imposed by B^T
    // (transpose) is not considered in the first-order conditions for minimum.
    error(true, "Method EquilibriumConstraints::Control::titrationOf is not supported yet.");
    data.controls.titrants.push_back(titrant);
    return *this;
}

auto EquilibriumConstraints::Control::titrationOfEither(String titrant1, String titrant2) -> Control&
{
    error(true, "Method EquilibriumConstraints::Control::titrationOfEither has not been implemented yet.");
    return *this;
}

//=================================================================================================
//
// EquilibriumConstraints::Until
//
//=================================================================================================

EquilibriumConstraints::Until::Until(EquilibriumConstraints::Data& data)
: data(data)
{}

auto EquilibriumConstraints::Until::volume(real value, String unit) -> Until&
{
    value = units::convert(value, unit, "m3");
    return custom([=](EquilibriumConstraintArgs args) { return args.props.volume() - value; });
}

auto EquilibriumConstraints::Until::internalEnergy(real value, String unit) -> Until&
{
    value = units::convert(value, unit, "J");
    return custom([=](EquilibriumConstraintArgs args) { return args.props.internalEnergy() - value; });
}

auto EquilibriumConstraints::Until::enthalpy(real value, String unit) -> Until&
{
    value = units::convert(value, unit, "J");
    return custom([=](EquilibriumConstraintArgs args) { return args.props.enthalpy() - value; });
}

auto EquilibriumConstraints::Until::gibbsEnergy(real value, String unit) -> Until&
{
    value = units::convert(value, unit, "J");
    return custom([=](EquilibriumConstraintArgs args) { return args.props.gibbsEnergy() - value; });
}

auto EquilibriumConstraints::Until::helmholtzEnergy(real value, String unit) -> Until&
{
    value = units::convert(value, unit, "J");
    return custom([=](EquilibriumConstraintArgs args) { return args.props.helmholtzEnergy() - value; });
}

auto EquilibriumConstraints::Until::entropy(real value, String unit) -> Until&
{
    value = units::convert(value, unit, "J/K");
    return custom([=](EquilibriumConstraintArgs args) { return args.props.entropy() - value; });
}

auto EquilibriumConstraints::Until::custom(const EquilibriumConstraintFn& fn) -> Until&
{
    data.fconstraints.push_back(fn);
    return *this;
}

//=================================================================================================
//
// EquilibriumConstraints::Preserve
//
//=================================================================================================

EquilibriumConstraints::Preserve::Preserve(EquilibriumConstraints::Data& data)
: data(data)
{}

auto EquilibriumConstraints::Preserve::volume() -> Preserve&
{
    return custom([](const ChemicalProps& props) { return props.volume(); });
}

auto EquilibriumConstraints::Preserve::internalEnergy() -> Preserve&
{
    return custom([](const ChemicalProps& props) { return props.internalEnergy(); });
}

auto EquilibriumConstraints::Preserve::enthalpy() -> Preserve&
{
    return custom([](const ChemicalProps& props) { return props.enthalpy(); });
}

auto EquilibriumConstraints::Preserve::gibbsEnergy() -> Preserve&
{
    return custom([](const ChemicalProps& props) { return props.gibbsEnergy(); });
}

auto EquilibriumConstraints::Preserve::helmholtzEnergy() -> Preserve&
{
    return custom([](const ChemicalProps& props) { return props.helmholtzEnergy(); });
}

auto EquilibriumConstraints::Preserve::entropy() -> Preserve&
{
    return custom([](const ChemicalProps& props) { return props.entropy(); });
}

auto EquilibriumConstraints::Preserve::custom(const ChemicalPropertyFn& fn) -> Preserve&
{
    error(!fn, "Imposing an empty custom chemical property function is not allowed.");
    data.pconstraints.push_back(fn);
    return *this;
}

//=================================================================================================
//
// EquilibriumConstraints::Fix
//
//=================================================================================================

EquilibriumConstraints::Fix::Fix(const ChemicalSystem& system, EquilibriumConstraints::Data& data)
: system(system), data(data)
{}

auto EquilibriumConstraints::Fix::chemicalPotential(const ChemicalFormula& substance, const Fn<real(real,real)>& fn) -> Fix&
{
    data.uconstraints.push_back({substance, fn});
    return *this;
}

auto EquilibriumConstraints::Fix::chemicalPotential(String substance, real value, String unit) -> Fix&
{
    value = units::convert(value, unit, "J/mol");
    return chemicalPotential(substance, [=](real T, real P) { return value; });
}

auto EquilibriumConstraints::Fix::lnActivity(const Species& species, real value) -> Fix&
{
    const auto R = universalGasConstant;

    auto fn = [=](real T, real P)
    {
        const auto u0 = species.props(T, P).G0;
        return u0 + R*T*value;
    };

    return chemicalPotential(species.formula(), fn);
}

auto EquilibriumConstraints::Fix::lnActivity(String name, real value) -> Fix&
{
    const auto idx = system.database().species().findWithName(name);
    error(idx >= system.database().species().size(),
        "Could not impose an activity constraint for species `", name, "` "
        "because it is not in the database.");
    const auto species = system.database().species()[idx];
    return lnActivity(species, value);
}

auto EquilibriumConstraints::Fix::lgActivity(String name, real value) -> Fix&
{
    return lnActivity(name, value * ln10);
}

auto EquilibriumConstraints::Fix::activity(String name, real value) -> Fix&
{
    return lnActivity(name, log(value));
}

auto EquilibriumConstraints::Fix::fugacity(String gas, real value, String unit) -> Fix&
{
    value = units::convert(value, unit, "bar");
    const auto& gases = system.database().speciesWithAggregateState(AggregateState::Gas);
    const auto idx = gases.findWithFormula(gas);
    error(idx >= gases.size(),
        "Could not impose a fugacity constraint for gas `", gas, "` because "
        "there is no gaseous species in the database with this chemical formula.");
    return lnActivity(gases[idx], value);
}

auto EquilibriumConstraints::Fix::pH(real value) -> Fix&
{
    const auto aqspecies = system.database().speciesWithAggregateState(AggregateState::Aqueous);
    const auto idx = aqspecies.findWithFormula("H+");
    error(idx >= aqspecies.size(),
        "Could not impose pH constraint because the database has "
        "no aqueous species with chemical formula `H+`.");
    return lnActivity(aqspecies[idx], -value * ln10); // pH = -log10(a[H+]) => ln(a[H+]) = -pH * ln10
}

auto EquilibriumConstraints::Fix::pMg(real value) -> Fix&
{
    const auto aqspecies = system.database().speciesWithAggregateState(AggregateState::Aqueous);
    const auto idx = aqspecies.findWithFormula("Mg+2");
    error(idx >= aqspecies.size(),
        "Could not impose pMg constraint because the database has "
        "no aqueous species with chemical formula `Mg+2`.");
    return lnActivity(aqspecies[idx], -value * ln10); // pMg = -log10(a[Mg+2]) => ln(a[Mg+2]) = -pH * ln10
}

auto EquilibriumConstraints::Fix::pe(real value) -> Fix&
{
    const auto ue0 = 0.0;            // the standard chemical potential of the electron species
    const auto lnae = -value * ln10; // the ln activity of the electron species
    const auto R = universalGasConstant;
    return chemicalPotential("e-", [=](real T, real P) { return ue0 + R*T*lnae; });
}

auto EquilibriumConstraints::Fix::Eh(real value, String unit) -> Fix&
{
    value = units::convert(value, unit, "V"); // in V = J/C
    const auto F = faradayConstant;           // in C/mol
    const auto ue = -F * value;               // in J/mol (chemical potential of electron)
    return chemicalPotential("e-", [=](real T, real P) { return ue; });
}

//=================================================================================================
//
// EquilibriumConstraints::Prevent
//
//=================================================================================================

EquilibriumConstraints::Prevent::Prevent(const ChemicalSystem& system, EquilibriumConstraints::Data& data)
: system(system), data(data)
{}

auto EquilibriumConstraints::Prevent::fromReacting(Index ispecies) -> Prevent&
{
    const auto size = system.species().size();
    error(ispecies >= size,
        "The given species index ", ispecies, " is out of bounds, "
        "since there are only ", size, " species in the chemical system.");
    data.rconstraints.species_cannot_react.insert(ispecies);
    return *this;
}

auto EquilibriumConstraints::Prevent::fromReacting(Pairs<Index, double> equation) -> Prevent&
{
    data.rconstraints.reactions_cannot_react.push_back(equation);
    return *this;
}

auto EquilibriumConstraints::Prevent::fromReacting(String what) -> Prevent&
{
    if(what.find("=") == String::npos) {
        return fromReacting(system.species().indexWithName(what));
    }
    else {
        Pairs<Index, double> pairs;
        for(auto [name, coeff] : parseReactionEquation(what))
            pairs.emplace_back(system.species().indexWithName(name), coeff);
        return fromReacting(pairs);
    }
}

auto EquilibriumConstraints::Prevent::fromIncreasing(Index ispecies) -> Prevent&
{
    const auto size = system.species().size();
    error(ispecies >= size,
        "The given species index ", ispecies, " is out of bounds, "
        "since there are only ", size, " species in the chemical system.");
    data.rconstraints.species_cannot_increase.insert(ispecies);
    return *this;
}

auto EquilibriumConstraints::Prevent::fromIncreasing(String species) -> Prevent&
{
    return fromIncreasing(system.species().indexWithName(species));
}

auto EquilibriumConstraints::Prevent::fromDecreasing(Index ispecies) -> Prevent&
{
    const auto size = system.species().size();
    error(ispecies >= size,
        "The given species index ", ispecies, " is out of bounds, "
        "since there are only ", size, " species in the chemical system.");
    data.rconstraints.species_cannot_decrease.insert(ispecies);
    return *this;
}

auto EquilibriumConstraints::Prevent::fromDecreasing(String species) -> Prevent&
{
    return fromDecreasing(system.species().indexWithName(species));
}

} // namespace Reaktoro
