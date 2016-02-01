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

class ChemicalPlot
{
public:
    ChemicalPlot();

    explicit ChemicalPlot(const ChemicalSystem& system);

    explicit ChemicalPlot(const ReactionSystem& reactions);

    virtual ~ChemicalPlot();

    auto name(std::string name) -> void;

    auto xdata(std::string x) -> void;

    auto ydata(std::vector<std::string> y) -> void;

    auto ydata(std::string y) -> void;

    auto xlabel(std::string) -> void;

    auto ylabel(std::string) -> void;

    auto xtics(std::string) -> void;

    auto ytics(std::string) -> void;

    auto xformat(std::string) -> void;

    auto yformat(std::string) -> void;

    auto xlogscale(int base=10) -> void;

    auto ylogscale(int base=10) -> void;

    auto legend(std::vector<std::string> legend) -> void;

    auto legend(std::string legend) -> void;

    auto nolegend() -> void;

    auto key(std::string) -> void;

    auto frequency(unsigned frequency) -> void;

    auto operator<<(std::string command) -> ChemicalPlot&;

    auto operator<<(std::stringstream command) -> ChemicalPlot&;

    auto open() -> void;

    auto update(const ChemicalState& state, double t) -> void;

    /// Compare a ChemicalPlot instance for equality
    auto operator==(const ChemicalPlot& other) -> bool;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};


} // namespace Reaktoro
