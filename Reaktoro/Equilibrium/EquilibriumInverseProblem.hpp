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
#include <map>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/ScalarTypes.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProperties;
class ChemicalState;
class ChemicalSystem;
class Partition;
class Phase;
class Species;
struct EquilibriumOptions;
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

/// A type used to represent a titrant and its attributes.
class Titrant
{
public:
    /// Construct a default Titrant instance.
    Titrant();

    /// Construct a Titrant instance from a Species instance.
    explicit Titrant(const Species& species);

    /// Construct a Titrant instance from a compound name.
    explicit Titrant(std::string compound);

    /// Construct a Titrant instance from a compound name.
    /// @param name The name of the titrant.
    /// @param formula The elemental formula of the titrant.
    Titrant(std::string name, const std::map<std::string, double>& formula);

    /// Return the name of the titrant.
    auto name() const -> std::string;

    /// Return the elemental formula of the titrant.
    auto formula() const -> std::map<std::string, double>;

    /// Return the elemental formula vector of the titrant.
    /// This method returns the vector of element stoichiometry of the titrant
    /// with respect to the elements in a chemical system.
    auto formula(const ChemicalSystem& system) const -> Vector;

    /// Return the molar mass of the titrant (in units of kg/mol).
    auto molarMass() const -> double;

private:
    /// The name of the titrant.
    std::string _name;

    /// The elemental formula of the titrant as pairs (element, stoichiometry).
    std::map<std::string, double> _formula;

    /// The molar mass of the titrant (in units of kg/mol)
    double _molar_mass;
};

/// Compare two Titrant instances for less than
auto operator<(const Titrant& lhs, const Titrant& rhs) -> bool;

/// Compare two Titrant instances for equality
auto operator==(const Titrant& lhs, const Titrant& rhs) -> bool;

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
    /// Construct an EquilibriumInverseProblem instance
    explicit EquilibriumInverseProblem(const ChemicalSystem& system);

    /// Construct a copy of an EquilibriumInverseProblem instance
    EquilibriumInverseProblem(const EquilibriumInverseProblem& other);

    /// Destroy this EquilibriumInverseProblem instance
    virtual ~EquilibriumInverseProblem();

    /// Assign a copy of an EquilibriumInverseProblem instance
    auto operator=(EquilibriumInverseProblem other) -> EquilibriumInverseProblem&;

    /// Add an activity constraint to the inverse equilibrium problem.
    /// @param species The name of the species for which its activity is given.
    /// @param value The value of the species activity.
    auto addSpeciesActivityConstraint(std::string species, double value) -> void;

    /// Add an amount constraint to the inverse equilibrium problem.
    /// @param species The name of the species for which its amount is given.
    /// @param value The value of the species amount (in units of mol).
    auto addSpeciesAmountConstraint(std::string species, double value) -> void;

    /// Add a phase amount constraint to the inverse equilibrium problem.
    /// @param phase The name of the phase for which its total amount is given.
    /// @param value The value of the phase total amount (in units of mol)
    auto addPhaseAmountConstraint(std::string phase, double value) -> void;

    /// Add a phase volume constraint to the inverse equilibrium problem.
    /// @param phase The name of the phase for which its volume is given.
    /// @param value The value of the phase volume (in units of m3)
    auto addPhaseVolumeConstraint(std::string phase, double value) -> void;

    /// Add a sum of phase volumes constraint to the inverse equilibrium problem.
    /// @param phases The names of the phases for which their sum of volumes is given.
    /// @param value The value of the volume sum of the phases (in units of m3)
    auto addSumPhaseVolumesConstraint(const std::vector<std::string>& phases, double value) -> void;

    /// Set the initial known molar amounts of the elements in the equilibrium partition.
    /// These are the amounts of the equilibrium elements before unknown amounts of titrants are added.
    auto setElementInitialAmounts(const Vector& b0) -> void;

    /// Set the initial amount of a titrant as an initial guess for the inverse equilibrium calculation.
    /// The titrant must have been added before using one of the `addTitrant` methods.
    /// @param titrant The name of the titrant.
    /// @param amount The amount of the titrant.
    /// @see addTitrant, addTitrants
    auto setTitrantInitialAmount(std::string titrant, double amount) -> void;

    /// Add a titrant to the inverse equilibrium problem.
    auto addTitrant(const Titrant& titrant) -> void;

    /// Add a titrant to the inverse equilibrium problem using either a species name or a compound.
    /// If `titrant` is the name of a species in the chemical system, then this species attributes
    /// are used to create a Titrant instance. If not, the string `titrant` is parsed to obtain its
    /// elemental composition and molar mass to create the needed Titrant instance.
    auto addTitrant(std::string titrant) -> void;

    /// Set two titrants to be mutually exclusive.
    /// Mutually exclusive titrants are used when only one of the titrants is allowed
    /// to be non-zero, while the other is zero. For example, one might specify the pH of
    /// the solution, whose value could be achieved by either the addition of an acid HCl
    /// or a base NaOH. To avoid infinitely many numerical solutions, it is important that
    /// these two titrants are mutually exclusive. If their amounts are `x1` and `x2`, an
    /// additional constraint is added so that `x1*x2 = 0`.
    auto setAsMutuallyExclusive(std::string titrant1, std::string titrant2) -> void;

    /// Return true if the instance has no equilibrium constraints.
    auto empty() const -> bool;

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
    auto residualEquilibriumConstraints(const Vector& x, const ChemicalState& state) const -> ResidualEquilibriumConstraints;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
