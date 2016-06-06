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
#include <istream>
#include <map>
#include <memory>
#include <string>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Kinetics/KineticState.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>

namespace Reaktoro {

/// A type used to describe the internal state of an interpreter.
struct InterpreterState
{
    /// The database used to initialize the chemical system
    Database database;

    /// The chemical editor used to define the chemical system
    ChemicalEditor editor;

    /// The chemical system used for the calculations
    ChemicalSystem system;

    /// The chemical reactions controlled by kinetics
    ReactionSystem reactions;

    /// The map of chemical states produced during the calculation
    std::map<std::string, KineticState> states;

    /// The defined mineral reactions
    std::vector<MineralReaction> mineral_reactions;

    /// The list of all compound names found in the script file.
    std::vector<std::string> compounds;

    /// The list of all elements that compose the compounds found in the script file.
    std::vector<std::string> elements;
};

/// A type used to define operations that interpret a Reaktoro script file.
class Interpreter
{
public:
    /// Construct a default Interpreter instance.
    Interpreter();

    /// Construct a copy of an Interpreter instance
    Interpreter(const Interpreter& other);

    /// Destroy this instance
    virtual ~Interpreter();

    /// Assign an Interpreter instance to this instance
    auto operator=(Interpreter other) -> Interpreter&;

    /// Execute a Reaktoro input script as string.
    auto execute(std::string str) -> void;

    /// Execute a Reaktoro input script as a file.
    auto execute(std::istream& stream) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
