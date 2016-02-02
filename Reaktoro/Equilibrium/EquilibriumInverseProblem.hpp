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
struct EquilibriumOptions;
struct EquilibriumResult;

/// A class used for defining an inverse equilibrium problem.
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
    auto addActivityConstraint(std::string species, double value) -> void;

    /// Add an amount constraint to the inverse equilibrium problem.
    /// @param species The name of the species for which its amount is given.
    /// @param value The value of the species amount (in units of mol).
    auto addAmountConstraint(std::string species, double value) -> void;

    /// Add a phase amount constraint to the inverse equilibrium problem.
    /// @param phase The name of the phase for which its total amount is given.
    /// @param value The value of the phase total amount (in units of mol)
    auto addPhaseAmountConstraint(std::string phase, double value) -> void;

    /// Add a phase volume constraint to the inverse equilibrium problem.
    /// @param phase The name of the phase for which its volume is given.
    /// @param value The value of the phase volume (in units of m3)
    auto addPhaseVolumeConstraint(std::string phase, double value) -> void;

    /// Set the known molar amounts of the elements in the equilibrium partition.
    /// These are the amounts of the equilibrium elements before unknown amounts
    /// of titrants are added.
    auto setInitialElementAmounts(const Vector& b0) -> void;

    /// Add a titrant to the inverse equilibrium problem.
    /// The amount of the titrant is unknown, but determined from given equilibrium constraints.
    /// The formula of the titrant must be specified in this function. For example, HCl and
    /// CO2 would be specified as:
    /// ~~~
    /// problem.addTitrant("HCl", {{"H", 1.0}, {"Cl", 1.0}});
    /// problem.addTitrant("CO2", {{"C", 1.0}, {"O", 2.0}});
    /// ~~~
    /// @param titrant The name of the titrant.
    /// @param formula The formula of the titrant.
    auto addTitrant(std::string titrant, std::map<std::string, double> formula) -> void;

    /// Set two titrants to be mutually exclusive.
    /// Two mutually exclusive titrants are needed when only one of them is allowed to be
    /// positive, while the other is zero. For example, one might specify the pH of the solution,
    /// whose value could be achieved by either the addition of an acid HCl or a base NaOH. To avoid
    /// infinitely many numerical solutions, it is important that these two titrants
    /// are mutually exclusive.
    auto setAsMutuallyExclusive(std::string titrant1, std::string titrant2) -> void;

    /// Return the chemical system.
    auto system() const -> const ChemicalSystem&;

    /// Return the formula matrix of the titrants.
    /// The formula matrix of the titrants is defined as the matrix whose (j,i)th entry
    /// contains the stoichiometric coefficient of jth element in the ith titrant.
    /// Its dimension is `E x T`, where `T` is the number of titrants.
    auto formulaMatrixTitrants() const -> Matrix;

    /// Return the initial amounts of elements in the equilibrium partition.
    /// These are the values of element amounts in the equilibrium partition
    /// before unknown amounts of titrants are added.
    auto initialElementAmounts() const -> Vector;

    /// Return the residual of the equilibrium constraints and their partial molar derivatives.
    /// @param properties The chemical properties of the system.
    auto residualConstraints(const ChemicalProperties& properties) const -> ChemicalVector;

    /// Return the residual of the mutually exclusive constraints and their x-derivatives.
    /// @param x The amounts of the titrants.
    auto residualMutuallyExclusiveConstraints(const Vector& x) const -> ChemicalVector;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
