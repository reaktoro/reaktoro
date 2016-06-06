// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class EquilibriumState;
class Partition;
struct EquilibriumResult;

/// A type used to define the result of the evaluation of a system of equilibrium constraints.
struct ResidualEquilibriumConstraints
{
    /// The residual values of the equilibrium constraints.
    Vector val;

    /// The partial derivatives of the residuals w.r.t. titrant amounts x.
    Matrix ddx;

    /// The partial derivatives of the residuals w.r.t. species amounts n.
    Matrix ddn;
};

/// A class used for defining an inverse equilibrium problem.
/// In an inverse equilibrium problem, not all elements have known molar amounts.
/// Their amount constraints are replaced by other equilibrium constraints such as
/// fixed species amount or activity, or the volume or total amount of a phase.
/// Since the amounts of elements are not known a priori, an inverse equilibrium
/// calculation tries to determine amounts of titrants that can control the specified
/// equilibrium constraints. The amount of the titrants are unknown, and its addition
/// or removal is done over the calculation so that the equilibrium state is driven
/// towards a state where all given equilibrium constraints are satisfied.
class EquilibriumInverseProblem
{
public:
    /// Construct an EquilibriumInverseProblem instance.
    explicit EquilibriumInverseProblem(const ChemicalSystem& system);

    /// Construct a copy of an EquilibriumInverseProblem instance.
    EquilibriumInverseProblem(const EquilibriumInverseProblem& other);

    /// Destroy this EquilibriumInverseProblem instance.
    virtual ~EquilibriumInverseProblem();

    /// Assign a copy of an EquilibriumInverseProblem instance.
    auto operator=(EquilibriumInverseProblem other) -> EquilibriumInverseProblem&;

    /// Set the partition of the chemical system.
    /// Use this method to specify the equilibrium, kinetic, and inert species.
    auto setPartition(const Partition& partition) -> EquilibriumInverseProblem&;

    /// Set the temperature for the equilibrium calculation (in units of K).
    /// By default, the temperature is 25 &deg;C.
    /// @param val The temperature value (in units of K).
    auto setTemperature(double val) -> EquilibriumInverseProblem&;

    /// Set the temperature for the equilibrium calculation with given units.
    /// By default, the temperature is 25 &deg;C.
    /// @param val The temperature value.
    /// @param units The units of the temperature (K, degC, degF, degR, kelvin, celsius, fahrenheit, rankine).
    auto setTemperature(double val, std::string units) -> EquilibriumInverseProblem&;

    /// Set the pressure for the equilibrium calculation (in units of Pa).
    /// By default, the pressure is 1 bar.
    /// @param val The pressure value (in units of Pa).
    /// @param units The units of the pressure (K, degC, degF, degR, kelvin, celsius, fahrenheit, rankine).
    auto setPressure(double val) -> EquilibriumInverseProblem&;

    /// Set the pressure for the equilibrium calculation (in units of Pa).
    /// By default, the pressure is 1 bar.
    /// @param val The pressure value.
    /// @param units The units of the pressure (Pa, kPa, MPa, GPa, atm, mmHg, inHg, psi, kpsi, Mpsi, psf, bar, torr, inH2O, ftH2O, pascal).
    auto setPressure(double val, std::string units) -> EquilibriumInverseProblem&;

    /// Set the initial known molar amounts of the elements in the equilibrium partition.
    /// These are the amounts of the equilibrium elements before unknown amounts of titrants are added.
    auto setElementInitialAmounts(const Vector& values) -> EquilibriumInverseProblem&;

    /// Add a given amount of a compound or species to the initial equilibrium recipe.
    /// @param name The name of the compound or species
    /// @param amount The amount of the compound or species
    /// @param units The units of the amount (must be convertible to either mol or kg)
    auto add(std::string name, double amount, std::string units) -> EquilibriumInverseProblem&;

    /// Fix the molar amount of a species at equilibrium.
    /// @param species The name of the species for which its amount is given.
    /// @param value The value of the species amount (in units of mol).
    auto fixSpeciesAmount(std::string species, double value, std::string units) -> EquilibriumInverseProblem&;

    /// Fix the molar amount of a species at equilibrium with given titrant.
    /// @param species The name of the species for which its amount is given.
    /// @param value The value of the species amount (in units of mol).
    /// @param titrant The titrant that controls the species amount.
    auto fixSpeciesAmount(std::string species, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&;

    /// Fix the mass of a species at equilibrium.
    /// @param species The name of the species for which its amount is given.
    /// @param value The value of the species mass (in units of kg).
    auto fixSpeciesMass(std::string species, double value, std::string units) -> EquilibriumInverseProblem&;

    /// Fix the mass of a species at equilibrium with given titrant.
    /// @param species The name of the species for which its amount is given.
    /// @param value The value of the species amount (in units of kg).
    /// @param titrant The titrant that controls the species mass.
    auto fixSpeciesMass(std::string species, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&;

    /// Fix the activity of a species at equilibrium.
    /// @param species The name of the species for which its activity is given.
    /// @param value The value of the species activity.
    auto fixSpeciesActivity(std::string species, double value) -> EquilibriumInverseProblem&;

    /// Fix the activity of a species at equilibrium with given titrant.
    /// @param species The name of the species for which its activity is given.
    /// @param value The value of the species activity.
    /// @param titrant The titrant that controls the species activity.
    auto fixSpeciesActivity(std::string species, double value, std::string titrant) -> EquilibriumInverseProblem&;

    /// Fix the activity of a species at equilibrium with either one of the two given titrants.
    /// @param species The name of the species for which its activity is given.
    /// @param value The value of the species activity.
    /// @param titrant1 The first titrant that might control the species activity.
    /// @param titrant2 The second titrant that might control the species activity.
    auto fixSpeciesActivity(std::string species, double value, std::string titrant1, std::string titrant2) -> EquilibriumInverseProblem&;

    /// Fix the fugacity of a gaseous species.
    /// @param species The name of the gaseous species.
    /// @param value The value of the species fugacity.
    /// @param units The units of the fugacity value.
    auto fixSpeciesFugacity(std::string species, double value, std::string units) -> EquilibriumInverseProblem&;

    /// Fix the fugacity of a gaseous species with given titrant.
    /// @param species The name of the gaseous species.
    /// @param value The value of the species fugacity.
    /// @param units The units of the fugacity value.
    /// @param titrant The titrant that controls the fugacity value.
    auto fixSpeciesFugacity(std::string species, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&;

    /// Fix the total molar amount of a phase at equilibrium with given titrant.
    /// @param phase The name of the phase for which its total molar amount is given.
    /// @param value The value of the phase total amount (in units of mol)
    /// @param titrant The titrant that controls the phase amount.
    auto fixPhaseAmount(std::string phase, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&;

    /// Fix the volume of a phase at equilibrium with given titrant.
    /// @param phase The name of the phase for which its volume is given.
    /// @param value The value of the phase volume (in units of m3)
    /// @param titrant The titrant that controls the phase volume.
    auto fixPhaseVolume(std::string phase, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&;

    /// Fix the total volume of a set of phases at equilibrium with given titrant.
    /// @param phases The names of the phases composing the phase set.
    /// @param value The value of the total volume of the phase set (in units of m3)
    /// @param titrant The titrant that controls the total volume of the phase set.
    auto fixPhaseSetVolume(const std::vector<std::string>& phases, double value, std::string units, std::string titrant) -> EquilibriumInverseProblem&;

    /// Fix the pH of the aqueous solution.
    /// @param value The pH value of the aqueous solution.
    auto pH(double value) -> EquilibriumInverseProblem&;

    /// Fix the pH of the aqueous solution with given titrant.
    /// @param value The pH value of the aqueous solution.
    /// @param titrant The titrant that control the solution pH.
    auto pH(double value, std::string titrant) -> EquilibriumInverseProblem&;

    /// Fix the pH of the aqueous solution with either one of the two given titrants.
    /// @param value The pH value of the aqueous solution.
    /// @param titrant1 The first titrant that might control the solution pH.
    /// @param titrant2 The second titrant that might control the solution pH.
    auto pH(double value, std::string titrant1, std::string titrant2) -> EquilibriumInverseProblem&;

    /// Fix the pe of the aqueous solution.
    /// @param value The pe value of the aqueous solution.
    auto pe(double value) -> EquilibriumInverseProblem&;

    /// Fix the pe of the aqueous solution with given half reaction.
    /// @param value The pe value of the aqueous solution.
    /// @param reaction The half reaction from which pe should be calculated.
    auto pe(double value, std::string reaction) -> EquilibriumInverseProblem&;

    /// Fix the Eh of the aqueous solution (in units of volts).
    /// @param value The Eh value of the aqueous solution.
    auto Eh(double value) -> EquilibriumInverseProblem&;

    /// Fix the Eh of the aqueous solution with given half reaction.
    /// @param value The Eh value of the aqueous solution.
    /// @param reaction The half reaction from which Eh should be calculated.
    auto Eh(double value, std::string reaction) -> EquilibriumInverseProblem&;

    /// Return a reference to the ChemicalSystem instance used to create this EquilibriumProblem instance
    auto system() const -> const ChemicalSystem&;

    /// Return a reference to the Partition instance used to create this EquilibriumProblem instance
    auto partition() const -> const Partition&;

    /// Return the temperature for the equilibrium calculation (in units of K)
    auto temperature() const -> double;

    /// Return the pressure for the equilibrium calculation (in units of Pa)
    auto pressure() const -> double;

    /// Return the number of constraints used in the inverse equilibrium problem.
    auto numConstraints() const -> Index;

    /// Return the number of titrants used in the inverse equilibrium problem.
    auto numTitrants() const -> Index;

    /// Return the formula matrix of the titrants.
    /// The formula matrix of the titrants is defined as the matrix whose (j,i)th entry
    /// contains the stoichiometric coefficient of jth element in the ith titrant.
    /// Its dimension is `E x T`, where `T` is the number of titrants.
    auto formulaMatrixTitrants() const -> Matrix;

    /// Return the initial amounts of elements in the equilibrium partition.
    /// These are the values of element amounts in the equilibrium partition
    /// before unknown amounts of titrants are added.
    auto elementInitialAmounts() const -> Vector;

    /// Return the initial amounts of titrants as initial guess for the inverse equilibrium calculation.
    auto titrantInitialAmounts() const -> Vector;

    /// Return the residuals of the equilibrium constraints and their partial derivatives.
    /// @param x The amounts of the titrants (in units of mol)
    /// @param state The chemical state of the system
    auto residualEquilibriumConstraints(const Vector& x, const EquilibriumState& state) const -> ResidualEquilibriumConstraints;

    /// Solve the inverse equilibrium problem.
    /// @param state The initial guess for the final chemical state solution.
    auto solve(EquilibriumState& state) -> EquilibriumResult;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
