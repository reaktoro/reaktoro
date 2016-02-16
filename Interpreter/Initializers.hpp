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
#include <string>

namespace Reaktoro {

// Forward declarations
class ChemicalEditor;
class ChemicalPlot;
class ChemicalSystem;
class EquilibriumPath;
class EquilibriumProblem;
class KineticPath;
class MineralReaction;
class ReactionSystem;

// Forward declarations of keyword types
namespace kwd {

class EquilibriumPath;
class EquilibriumProblem;
class KineticPath;
class MineralReaction;
class Plot;

} // namespace kwd

/// Initialize a ChemicalEditor object using the whole input script file.
/// This method is used to initialize a ChemicalEditor instance by firstly
/// identifying all compound and species names in the input file. After this,
/// an aqueous phase is created with all possible species found in the database
/// that contains the elements composing the list of found compounds. A gaseous
/// phase is created only if names of gaseous species present in the database are
/// found in the input script. Pure mineral phases are created for each mineral name
/// found in the input script file that is also present in the database.
auto initializeChemicalEditor(ChemicalEditor& editor, std::string script) -> void;

/// Initialize a ChemicalPlot object using a kwd::Plot object.
auto initializeChemicalPlot(ChemicalPlot& plot, const kwd::Plot& keyword) -> void;

/// Initialize a MineralReaction object using a kwd::MineralReaction object.
auto initializeMineralReaction(MineralReaction& reaction, const kwd::MineralReaction& keyword) -> void;

/// Initialize an EquilibriumProblem object using a kwd::EquilibriumProblem object.
auto initializeEquilibriumProblem(EquilibriumProblem& problem, const kwd::EquilibriumProblem& keyword) -> void;

/// Initialize an EquilibriumPath object using a kwd::EquilibriumPath object.
auto initializeEquilibriumPath(EquilibriumPath& path, const kwd::EquilibriumPath& keyword) -> void;

/// Initialize a KineticPath object using a kwd::KineticPath object.
auto initializeKineticPath(KineticPath& path, const kwd::KineticPath& keyword) -> void;

} // namespace Reaktoro
