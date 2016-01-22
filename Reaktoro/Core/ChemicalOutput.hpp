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

class ChemicalOutput
{
public:
    ChemicalOutput();

    explicit ChemicalOutput(const ChemicalSystem& system);

    explicit ChemicalOutput(const ReactionSystem& reactions);

    virtual ~ChemicalOutput();

    auto file(std::string filename) -> void;

    auto terminal(bool active) -> void;

    auto data(std::vector<std::string> quantities) -> void;

    auto data(std::string quantities) -> void;

    auto header(std::vector<std::string> header) -> void;

    auto header(std::string header) -> void;

    auto open() -> void;

    auto update(const ChemicalState& state, double t) -> void;

    auto close() -> void;

    operator bool() const;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
