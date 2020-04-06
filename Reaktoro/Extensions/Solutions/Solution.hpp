// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#pragma once

// C++ includes
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// A type used to describe the state of a solution.
struct SolutionState
{
    /// The temperature of the solution (in K).
    real T = {};

    /// The pressure of the solution (in Pa).
    real P = {};

    /// The mole fractions of the species in the solution.
    VectorXr x;
};

/// Compare two SolutionState instances for equality.
auto operator==(const SolutionState& l, const SolutionState& r) -> bool;

/// The base class for all derived solution classes.
/// @ingroup Solutions
class Solution
{
public:
    /// Construct a default Solution instance.
    Solution();

    /// Construct a Solution instance with given species.
    /// @param species The names of the species in the solution
    explicit Solution(std::vector<Species> species);

    /// Set the name of the solution.
    auto setName(std::string name) -> void;

    /// Return the number of species in the solution
    auto numSpecies() const -> unsigned;

    /// Return the name of the solution.
    auto name() const -> std::string;

    /// Return the species that compose the solution
    /// @return The species that compose the solution
    auto species() const -> const std::vector<Species>&;

    /// Return a species in the solution
    /// @param index The index of the species
    /// @return The species with given index
    auto species(const Index& index) const -> const Species&;

    /// Return the index of a species in the solution
    /// @param name The name of the species in the solution
    /// @return The index of the species if found, or the number of species otherwise
    auto indexSpecies(const std::string& name) const -> Index;

    /// Return the index of the first species in the solution with any of the given names.
    /// @param names The tentative names of the species in the solution.
    /// @return The index of the species if found, or the number of species otherwise.
    auto indexSpeciesAny(const std::vector<std::string>& names) const -> Index;

    /// Return the names of the species in the solution
    auto namesSpecies() const -> std::vector<std::string>;

    /// Return the charges of the species in the solution
    auto charges() const -> VectorXr;

    /// Calculates the mole fractions of the species and their partial derivatives
    /// @param n The molar abundance of the species (in mol)
    /// @return The mole fractions and their partial derivatives
    auto moleFractions(VectorXrConstRef n) const -> VectorXr;

    /// Calculate the state of the solution.
    /// @param T The temperature (in K)
    /// @param P The pressure (in Pa)
    /// @param n The molar amounts of the species in the solution (in mol)
    auto state(real T, real P, VectorXrConstRef n) const -> SolutionState;

private:
    /// The name of solution
    std::string m_name;

    /// The species in the solution
    std::vector<Species> m_species;
};

} // namespace Reaktoro
