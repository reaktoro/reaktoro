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
#include <vector>
#include <sstream>
#include <string>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class ReactionSystem;

/// A class used to create plots from sequence of chemical states.
class ChemicalPlot
{
public:
    /// Construct a default ChemicalPlot instance.
    ChemicalPlot();

    /// Construct a ChemicalPlot instance using a ChemicalSystem instance.
    explicit ChemicalPlot(const ChemicalSystem& system);

    /// Construct a ChemicalPlot instance using a ReactionSystem instance.
    explicit ChemicalPlot(const ReactionSystem& reactions);

    /// Destroy this ChemicalPlot instance.
    virtual ~ChemicalPlot();

    /// Set the name of the plot.
    auto setName(std::string name) -> void;

    /// Set the chemical quantity that is plotted along the x-axis.
    auto setXData(std::string quantity) -> void;

    /// Add a chemical quantity to be plotted along the y-axis.
    auto addYData(std::string quantity) -> void;

    /// Add a chemical quantity to be plotted along the y-axis.
    auto addYData(std::string quantity, std::string legend) -> void;

    /// Set the label of the x-axis.
    auto setXLabel(std::string) -> void;

    /// Set the label of the y-axis.
    auto setYLabel(std::string) -> void;

    /// Set the tics of the x-axis.
    auto setXTics(std::string) -> void;

    /// Set the tics of the y-axis.
    auto setYTics(std::string) -> void;

    /// Set the numeric display format of the x-axis.
    auto setXFormat(std::string) -> void;

    /// Set the numeric display format of the y-axis.
    auto setYFormat(std::string) -> void;

    /// Set the x-axis to log-scale.
    auto setXLogscale(int base=10) -> void;

    /// Set the y-axis to log-scale.
    auto setYLogscale(int base=10) -> void;

    /// Set the key options.
    auto setKey(std::string) -> void;

    /// Set the refresh rate of the real-time plot.
    auto setRefreshRate(unsigned frequency) -> void;

    /// Enable or disable the addLegend.
    auto enableLegend(bool enable) -> void;

    /// Inject a gnuplot command to the script file.
    auto operator<<(std::string command) -> ChemicalPlot&;

    /// Inject a gnuplot command to the script file.
    auto operator<<(std::stringstream command) -> ChemicalPlot&;

    /// Open the plot.
    auto open() -> void;

    /// Update the plot with a new chemical state and a tag.
    auto update(const ChemicalState& state, double t) -> void;

    /// Compare a ChemicalPlot instance for equality
    auto operator==(const ChemicalPlot& other) -> bool;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
