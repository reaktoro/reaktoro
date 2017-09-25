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
class StringList;

/// A type used to output sequence of chemical states to a file or terminal.
class ChemicalOutput
{
public:
    /// Construct a default ChemicalOutput instance.
    ChemicalOutput();

    /// Construct a ChemicalOutput instance with given ChemicalSystem instance.
    explicit ChemicalOutput(const ChemicalSystem& system);

    /// Construct a ChemicalOutput instance with given ReactionSystem instance.
    explicit ChemicalOutput(const ReactionSystem& reactions);

    /// Destroy this ChemicalOutput instance.
    virtual ~ChemicalOutput();

    /// Set the name of the output file.
    auto filename(std::string filename) -> void;

    /// Return the name of the output file.
    auto filename() const -> std::string;

    /// Add a quantity to be output.
    /// @param quantity The quantity name.
    auto add(std::string quantity) -> void;

    /// Add a quantity to be output.
    /// @param quantity The quantity name.
    /// @param label The label to be used in the headings.
    auto add(std::string quantity, std::string label) -> void;

    /// Set the titles of extra columns in the output file.
    auto attachments(std::vector<std::string> titles) -> void;

    /// Attach value to extra column data in the output file.
    auto attach(int value) -> void;

    /// Attach extra column data in the output file.
    auto attach(double value) -> void;

    /// Attach extra column data in the output file.
    auto attach(std::string value) -> void;

    /// Set the floating-point precision in the output.
    auto precision(int val) -> void;

    /// Enable or disable output in scientific format.
    auto scientific(bool enable) -> void;

    /// Enable or disable the output to the terminal.
    auto terminal(bool enabled) -> void;

    /// Return the name of the quantities in the output file.
    auto quantities() const -> std::vector<std::string>;

    /// Return the headings of the output file.
    auto headings() const -> std::vector<std::string>;

    /// Open the output file.
    auto open() -> void;

    /// Update the output with a new chemical state and its tag.
    auto update(const ChemicalState& state, double t) -> void;

    /// Close the output file.
    auto close() -> void;

    /// Convert this ChemicalOutput instance to bool.
    operator bool() const;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
