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
#include <vector>
#include <sstream>
#include <string>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;
class ReactionSystem;

struct ChemicalPlotOptions
{
    /// The name of the plot.
    std::string name;

    /// The quantity that spans the x-axis
    std::string x;

    /// The quantities to be plotted along the y-axis
    std::vector<std::string> y;

    /// The names of each curve given by member @ref y.
    std::vector<std::string> legend;

    /// The frequency in which the plot is refreshed per second.
    unsigned frequency = 30;

    /// The Gnuplot commands used to configure the plot.
    std::string config;
};

class ChemicalPlot
{
public:
    ChemicalPlot();

    explicit ChemicalPlot(const ChemicalSystem& system);

    explicit ChemicalPlot(const ReactionSystem& reactions);

    virtual ~ChemicalPlot();

    auto name(std::string name) -> void;

    auto x(std::string x) -> void;

    auto y(std::vector<std::string> y) -> void;

    auto y(std::string y) -> void;

    auto legend(std::vector<std::string> legend) -> void;

    auto legend(std::string legend) -> void;

    auto frequency(unsigned frequency) -> void;

    auto operator<<(std::string command) -> ChemicalPlot&;

    auto operator<<(std::stringstream command) -> ChemicalPlot&;

    auto open() -> void;

    auto update(const ChemicalState& state, double t) -> void;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
