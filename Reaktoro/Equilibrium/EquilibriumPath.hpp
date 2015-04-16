// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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
#include <string>

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class ChemicalState;
class Partition;

/// A struct that describes the options from an equilibrium path calculation.
struct EquilibriumPathOptions
{
    /// The number of points that discretize the equilibrium path.
    unsigned num_points = 20;

    /// The options for the equilibrium calculations.
    EquilibriumOptions equilibrium;

    /// The string containing Gnuplot commands for customized plotting
    std::string gnuplot;
};

/// A class that describes a path of equilibrium states.
class EquilibriumPath
{
public:
    /// Construct a default EquilibriumPath instance
    EquilibriumPath();

    /// Construct an EquilibriumPath instance
    explicit EquilibriumPath(const ChemicalSystem& system);

    /// Construct a copy of an EquilibriumPath instance
    EquilibriumPath(const EquilibriumPath& other);

    /// Destroy an EquilibriumPath instance
    virtual ~EquilibriumPath();

    /// Assign an EquilibriumPath instance to this instance
    auto operator=(EquilibriumPath other) -> EquilibriumPath&;

    /// Set the options for the equilibrium path calculation and visualization
    auto setOptions(const EquilibriumPathOptions& options) -> void;

    /// Set the partition of the chemical system
    auto setPartition(const Partition& partition) -> void;

    /// Set the partition of the chemical system using a formatted string
    auto setPartition(std::string partition) -> void;

    /// Solve the path of equilibrium states between two chemical states
    auto solve(const ChemicalState& state_i, const ChemicalState& state_f) -> void;

    /// Output the equilibrium path
    auto output(std::string list) -> void;

    /// Plot the equilibrium path
    auto plot(std::string list) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
