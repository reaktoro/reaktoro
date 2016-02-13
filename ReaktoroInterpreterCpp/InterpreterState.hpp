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
#include <map>
#include <string>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
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
    std::map<std::string, ChemicalState> states;

    /// The defined mineral reactions
    std::vector<MineralReaction> mineral_reactions;
};

} // namespace Reaktoro
