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
#include <string>
#include <vector>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;
class Phreeqc;

class PhreeqcEditor
{
public:
	/// Construct a default PhreeqcEditor instance.
	PhreeqcEditor();

	/// Construct a PhreeqcEditor instance with given database file.
	/// @param database The path, including name, to the PHREEQC database file.
	PhreeqcEditor(std::string database);

	/// Construct a copy of a PhreeqcEditor instance.
	PhreeqcEditor(const PhreeqcEditor& other);

	/// Destroy this PhreeqcEditor instance
	virtual ~PhreeqcEditor();

	/// Assign another PhreeqcEditor instance to this.
	auto operator=(PhreeqcEditor other) -> PhreeqcEditor&;

	/// Set the PHREEQC database file to be used by the PhreeqcEditor.
	/// @param database The path, including name, to the PHREEQC database file.
	auto setDatabase(std::string database) -> void;

	/// Set the aqueous species in the system by specifying which elements should exist.
	/// @param elements The names of the elements.
	auto setAqueousPhase(const std::vector<std::string>& elements) -> void;

	/// Set the aqueous species in the system by specifying which elements should exist.
	/// @param elements The names of the elements as a space-separated string list.
	auto setAqueousPhase(std::string elements) -> void;

	/// Set the gaseous phase in the system by specifying the end-member gases.
	/// @param gases The names of the gases.
	auto setGaseousPhase(const std::vector<std::string>& gases) -> void;

	/// Set the gaseous phase in the system by specifying the end-member gases.
	/// @param gases The names of the gases as a space-separated string list.
	auto setGaseousPhase(std::string gases) -> void;

	/// Set the mineral phases in the system by specifying the names of the pure minerals.
	/// @param minerals The names of the pure minerals.
	auto setMineralPhases(const std::vector<std::string>& minerals) -> void;

	/// Set the mineral phases in the system by specifying the names of the pure minerals.
	/// @param minerals The names of the pure minerals.
	auto setMineralPhases(std::string minerals) -> void;

	/// Convert this PhreeqcEditor instance into a ChemicalSystem instance
	operator ChemicalSystem() const;

	/// Convert this PhreeqcEditor instance into a Phreeqc instance
	operator Phreeqc() const;

private:
	struct Impl;

	std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
