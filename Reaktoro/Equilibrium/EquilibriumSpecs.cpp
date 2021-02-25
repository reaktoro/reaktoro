// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionEquation.hpp>

namespace Reaktoro {
namespace {

/// Return a Species object in a Database with given formula and aggregate state
auto getSpecies(const Database& db, const String& formula, AggregateState aggstate) -> Species
{
    const auto selected = db.speciesWithAggregateState(aggstate);
    const auto idx = selected.findWithFormula(formula);
    if(idx < selected.size()) return selected[idx];
    else return Species();
}

/// Return a Species object in a Database with given formula and aqueous aggregate state.
auto getAqueousSpecies(const Database& db, const String& formula) -> Species
{
    return getSpecies(db, formula, AggregateState::Aqueous);
}

/// Return a Species object in a Database with given formula and gaseous aggregate state.
auto getGaseousSpecies(const Database& db, const String& formula) -> Species
{
    return getSpecies(db, formula, AggregateState::Gas);
}

} // namespace

EquilibriumSpecs::EquilibriumSpecs(const ChemicalSystem& system)
: m_system(system)
{
}

//=================================================================================================
//
// METHODS TO SPECIFY THERMODYNAMIC CONSTRAINTS
//
//=================================================================================================

auto EquilibriumSpecs::temperature() -> void
{
    unknownT = false;
    addParameter("T");
}

auto EquilibriumSpecs::pressure() -> void
{
    unknownP = false;
    addParameter("P");
}

auto EquilibriumSpecs::volume() -> void
{
    EquilibriumConstraintEquation constraint;
    constraint.name = "volume";
    auto V = addParameter("V");
    constraint.fn = [=](const ChemicalProps& props) { return props.volume() - V; };
    addConstraint(constraint);
}

auto EquilibriumSpecs::internalEnergy() -> void
{
    EquilibriumConstraintEquation constraint;
    constraint.name = "internalEnergy";
    auto U = addParameter("U");
    constraint.fn = [=](const ChemicalProps& props) { return props.internalEnergy() - U; };
    addConstraint(constraint);
}

auto EquilibriumSpecs::enthalpy() -> void
{
    EquilibriumConstraintEquation constraint;
    constraint.name = "enthalpy";
    auto H = addParameter("H");
    constraint.fn = [=](const ChemicalProps& props) { return props.enthalpy() - H; };
    addConstraint(constraint);
}

auto EquilibriumSpecs::gibbsEnergy() -> void
{
    EquilibriumConstraintEquation constraint;
    constraint.name = "gibbsEnergy";
    auto G = addParameter("G");
    constraint.fn = [=](const ChemicalProps& props) { return props.gibbsEnergy() - G; };
    addConstraint(constraint);
}

auto EquilibriumSpecs::helmholtzEnergy() -> void
{
    EquilibriumConstraintEquation constraint;
    constraint.name = "helmholtzEnergy";
    auto A = addParameter("A");
    constraint.fn = [=](const ChemicalProps& props) { return props.helmholtzEnergy() - A; };
    addConstraint(constraint);
}

auto EquilibriumSpecs::entropy() -> void
{
    EquilibriumConstraintEquation constraint;
    constraint.name = "entropy";
    auto S = addParameter("S");
    constraint.fn = [=](const ChemicalProps& props) { return props.entropy() - S; };
    addConstraint(constraint);
}

//=================================================================================================
//
// METHODS TO SPECIFY CHEMICAL POTENTIAL CONSTRAINTS
//
//=================================================================================================

auto EquilibriumSpecs::chemicalPotential(String substance) -> void
{
    const auto pid = "u[" + substance + "]";
    EquilibriumConstraintChemicalPotential constraint;
    constraint.name = pid;
    constraint.substance = substance;
    auto u = addParameter(pid);
    constraint.fn = [=](const ChemicalProps& props) { return u; };
    addConstraint(constraint);
}

/// Create a chemical potential constraint from activity constraint.
/// @param species The Species object in the database for which activity value is given.
/// @param param The input parameter used to store the given activity value.
/// @param convert_to_lna The function that converts the activity value to natural log (e.g., convert lg(a[Ca+2]) to ln(a[Ca+2]), or pH to ln(a[H+]))
/// @return EquilibriumConstraintChemicalPotential
 auto createActivityConstraint(const Species& species, Param param, const Fn<real(real)>& convert_to_lna) -> EquilibriumConstraintChemicalPotential
{
    EquilibriumConstraintChemicalPotential constraint;
    constraint.name = param.id();
    constraint.substance = species.formula();
    constraint.fn = [=](const ChemicalProps& props)
    {
        const auto T = props.temperature();
        const auto P = props.pressure();
        const auto R = universalGasConstant;
        const auto u0 = species.props(T, P).G0;
        const auto lna = convert_to_lna(param);
        return u0 + R*T*lna;
    };
    return constraint;
}

auto EquilibriumSpecs::lnActivity(const Species& species) -> void
{
    const auto pid = "lnActivity[" + species.name() + "]";
    const auto convert_to_lna = [](real paramvalue) { return paramvalue; };
    auto lna = addParameter(pid);
    const auto constraint = createActivityConstraint(species, lna, convert_to_lna);
    addConstraint(constraint);
}

auto EquilibriumSpecs::lnActivity(String name) -> void
{
    const auto specieslist = m_system.database().species();
    const auto idx = specieslist.findWithName(name);
    error(idx >= specieslist.size(),
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
    error(species.name().empty(),
        "Could not impose the fugacity constraint for gas `", gas, "` because "
        "there is no gaseous species in the database with this chemical formula.");
    const auto pid = "f[" + gas + "]";
    auto fgas = addParameter(pid);
    const auto convert_to_lna = [](real fgas) { return log(fgas); };
    const auto constraint = createActivityConstraint(species, fgas, convert_to_lna);
    addConstraint(constraint);
}

auto EquilibriumSpecs::pH() -> void
{
    const auto species = getAqueousSpecies(m_system.database(), "H+");
    error(species.name().empty(),
        "Could not impose pH constraint because the database has "
        "no aqueous species with chemical formula `H+`.");
    auto pH = addParameter("pH");
    const auto convert_to_lna = [](real pH) { return -pH * ln10; };
    const auto constraint = createActivityConstraint(species, pH, convert_to_lna);
    addConstraint(constraint);
}

auto EquilibriumSpecs::pMg() -> void
{
    const auto species = getAqueousSpecies(m_system.database(), "Mg+2");
    error(species.name().empty(),
        "Could not impose pMg constraint because the database has "
        "no aqueous species with chemical formula `Mg+2`.");
    auto pMg = addParameter("pMg");
    const auto convert_to_lna = [](real pMg) { return -pMg * ln10; };
    const auto constraint = createActivityConstraint(species, pMg, convert_to_lna);
    addConstraint(constraint);
}

auto EquilibriumSpecs::pE() -> void
{
    const auto ue0 = 0.0; // the standard chemical potential of the electron species assumed zero
    const auto R = universalGasConstant;
    EquilibriumConstraintChemicalPotential constraint;
    constraint.name = "pE";
    constraint.substance = ChemicalFormula("e-");
    auto pE = addParameter("pE");
    constraint.fn = [=](const ChemicalProps& props)
    {
        const auto T = props.temperature();
        const auto lnae = -pE * ln10; // the ln activity of the electron species
        return ue0 + R*T*lnae;
    };
    addConstraint(constraint);
}

auto EquilibriumSpecs::Eh() -> void
{
    const auto F = faradayConstant; // in C/mol
    EquilibriumConstraintChemicalPotential constraint;
    constraint.name = "Eh";
    constraint.substance = ChemicalFormula("e-");
    auto Eh = addParameter("Eh");
    constraint.fn = [=](const ChemicalProps& props)
    {
        return -F * Eh; // in J/mol (chemical potential of electron)
    };
    addConstraint(constraint);
}

//=================================================================================================
//
// METHODS TO SPECIFY HOW THE CHEMICAL SYSTEM IS OPEN
//
//=================================================================================================

auto EquilibriumSpecs::openTo(const ChemicalFormula& substance) -> void
{
    throwErrorIfTitrantHasBeenRegistered(substance);
    titrants_explicit.push_back(substance);
}

//=================================================================================================
//
// METHODS TO GET THE NUMBER OF INTRODUCED CONSTRAINTS, PARAMETERS, AND CONTROL VARIABLES
//
//=================================================================================================

auto EquilibriumSpecs::numParameters() const -> Index
{
    return m_params.size();
}

auto EquilibriumSpecs::numControlVariables() const -> Index
{
    return unknownT + unknownP + numTitrants();
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

auto EquilibriumSpecs::numConstraints() const -> Index
{
    return econstraints.size() + uconstraints.size();
}

auto EquilibriumSpecs::numConstraintsEquationType() const -> Index
{
    return econstraints.size();
}

auto EquilibriumSpecs::numConstraintsChemicalPotentialType() const -> Index
{
    return uconstraints.size();
}

//=================================================================================================
//
// METHODS TO GET THE NAMES OF INTRODUCED CONSTRAINTS, PARAMETERS, AND CONTROL VARIABLES
//
//=================================================================================================

auto EquilibriumSpecs::namesParameters() const -> Strings
{
    return vectorize(m_params, RKT_LAMBDA(x, x.id()));
}

auto EquilibriumSpecs::namesControlVariables() const -> Strings
{
    Strings names;
    if(unknownT) names.push_back("T");
    if(unknownP) names.push_back("P");
    return concatenate(names, namesTitrants());
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
    return concatenate(namesConstraintsEquationType(), namesConstraintsChemicalPotentialType());
}

auto EquilibriumSpecs::namesConstraintsEquationType() const -> Strings
{
    return vectorize(econstraints, RKT_LAMBDA(x, x.name));
}

auto EquilibriumSpecs::namesConstraintsChemicalPotentialType() const -> Strings
{
    return vectorize(uconstraints, RKT_LAMBDA(x, x.name));
}

//=================================================================================================
//
// METHODS TO ADD CONSTRAINTS AND PARAMETERS
//
//=================================================================================================

auto EquilibriumSpecs::addConstraint(const EquilibriumConstraintEquation& constraint) -> void
{
    const auto constraint_has_same_name = containsfn(econstraints, RKT_LAMBDA(x, x.name == constraint.name));
    error(constraint_has_same_name, "Cannot impose a new equation constraint with same name (", constraint.name, ").");
    error(constraint.name.empty(), "An equation constraint cannot have an empty name.");
    error(!constraint.fn, "Imposing an empty function for an equation constraint is not allowed.");
    econstraints.push_back(constraint);
}

auto EquilibriumSpecs::addConstraint(const EquilibriumConstraintChemicalPotential& constraint) -> void
{
    const auto constraint_has_same_name = containsfn(uconstraints, RKT_LAMBDA(x, x.name == constraint.name));
    error(constraint_has_same_name, "Cannot impose a new chemical potential constraint with same name (", constraint.name, ").");
    error(constraint.name.empty(), "A chemical potential constraint cannot have an empty name.");
    error(!constraint.fn, "Imposing an empty function for a chemical potential constraint is not allowed.");
    throwErrorIfTitrantHasBeenRegistered(constraint.substance);
    uconstraints.push_back(constraint);
    titrants_implicit.push_back(constraint.substance);
}

auto EquilibriumSpecs::addParameter(String param) -> Param&
{
    errorif(m_params.exists(param), "Could not add parameter with identifier `", param, "` to the "
        "list of input parameters for the equilibrium problem because "
        "another parameter already exists with same identifier.");
    return m_params.append(param, 0.0);
}

auto EquilibriumSpecs::addParameter(Param param) -> Param&
{
    errorif(m_params.exists(param.id()), "Could not add parameter with identifier `", param.id(), "` to the "
        "list of input parameters for the equilibrium problem because "
        "another parameter already exists with same identifier.");
    return m_params.append(param);
}

//=================================================================================================
//
// MISCELLANEOUS METHODS
//
//=================================================================================================

auto EquilibriumSpecs::system() const -> const ChemicalSystem&
{
    return m_system;
}

auto EquilibriumSpecs::params() const -> Params
{
    return m_params;
}

auto EquilibriumSpecs::isTemperatureUnknown() const -> bool
{
    return unknownT;
}

auto EquilibriumSpecs::isPressureUnknown() const -> bool
{
    return unknownP;
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

auto EquilibriumSpecs::constraintsEquationType() const -> Vec<EquilibriumConstraintEquation> const&
{
    return econstraints;
}

auto EquilibriumSpecs::constraintsChemicalPotentialType() const -> Vec<EquilibriumConstraintChemicalPotential> const&
{
    return uconstraints;
}

//=================================================================================================
//
// PRIVATE METHODS
//
//=================================================================================================

auto EquilibriumSpecs::throwErrorIfTitrantHasBeenRegistered(const ChemicalFormula& substance) const -> void
{
    const auto already_an_explicit_titrant = containsfn(titrants_explicit, RKT_LAMBDA(x, x.equivalent(substance)));
    errorif(already_an_explicit_titrant, "Cannot specify again that the chemical system is open to substance ", substance.str(), ". "
        "You have already explicitly specified that the chemical system is open to this substance.");
    const auto idx = indexfn(titrants_implicit, RKT_LAMBDA(x, x.equivalent(substance)));
    const auto already_an_implicit_titrant = idx < titrants_implicit.size();
    errorif(already_an_implicit_titrant, "Cannot specify again that the chemical system is open to substance ", substance.str(), ". "
        "You have implicitly specified that the chemical system is open to this substance "
        "when imposing the ", uconstraints[idx].name, " constraint.");
}

} // namespace Reaktoro
