// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <map>
#include <memory>
#include <string>
#include <fstream>

// Reaktoro includes
#include <Reaktoro/Common/Json.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;

/// Used to interpret json files containing defined calculations.
class Interpreter
{
public:
    /// Construct a default Interpreter instance.
    Interpreter();

    /// Construct a copy of an Interpreter instance.
    Interpreter(const Interpreter& other);

    /// Destroy this Interpreter instance.
    virtual ~Interpreter();

    /// Assing another Interpreter instance to this.
    auto operator=(Interpreter other) -> Interpreter&;

    /// Execute an input script.
    /// @param input The input json object.
    auto executeJsonObject(json input) -> void;

    /// Execute an input script.
    /// @param input The json-formatted input string.
    auto executeJsonString(std::string input) -> void;

    /// Execute an input script.
    /// @param input The name of the json-formatted input file.
    auto executeJsonFile(std::string input) -> void;

    /// Return the constructed chemical system.
    auto system() -> const ChemicalSystem&;

    /// Return all saved chemical states during the execution.
    auto states() -> const std::map<std::string, ChemicalState>&;

    /// Return the saved chemical state with given reference name.
    auto state(std::string reference) -> const ChemicalState&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
