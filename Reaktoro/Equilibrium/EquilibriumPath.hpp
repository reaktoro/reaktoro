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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class ChemicalState;
class Partition;

/// A struct that describes the options for outputting in an equilibrium path calculation.
struct EquilibriumPathOutputOptions
{
    /// The flag that indicates if output is active.
    bool active = false;

    /// The flag that indicates if output should be done at the terminal.
    bool terminal = true;

    /// The name of the file to which the output should be written.
    std::string file;

    /// The names of the quantities to be output.
    std::vector<std::string> data;

    /// The names of the quantities to appear as column header in the output.
    std::vector<std::string> header;
};

/// A struct that describes the options for plotting in an equilibrium path calculation.
struct EquilibriumPathPlotOptions
{
    /// The quantity that spans the x-axis
    std::string x;

    /// The names of the quantities that spans the x-axis
    std::vector<std::string> y;

    /// The Gnuplot commands used to configure the plot.
    std::string config;

    /// The frequency in which the plot is refreshed per second.
    unsigned frequency = 30;
};

/// A struct that describes the options from an equilibrium path calculation.
struct EquilibriumPathOptions
{
    /// The options for the chemical equilibrium calculations.
    EquilibriumOptions equilibrium;

    /// The options for outputting the equilibrium path calculation.
    EquilibriumPathOutputOptions output;

    /// The options for plotting the equilibrium path calculation.
    std::vector<ChemicalPlotOptions> plots;
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
