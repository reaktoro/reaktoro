// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright © 2014-2021 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #pragma once

// // C++ includes
// #include <memory>
// #include <string>
// #include <vector>

// namespace Reaktoro {

// // Forward declarations
// class ChemicalSystem;
// class Phreeqc;
// class StringList;

// class PhreeqcEditor
// {
// public:
// /// Construct a default PhreeqcEditor instance.
// PhreeqcEditor();

// /// Construct a PhreeqcEditor instance with given database file.
// /// @param database The path, including name, to the PHREEQC database file.
// PhreeqcEditor(std::string database);

// /// Construct a copy of a PhreeqcEditor instance.
// PhreeqcEditor(const PhreeqcEditor& other);

// /// Destroy this PhreeqcEditor instance
// virtual ~PhreeqcEditor();

// /// Assign another PhreeqcEditor instance to this.
// auto operator=(PhreeqcEditor other) -> PhreeqcEditor&;

// /// Set the PHREEQC database file to be used by the PhreeqcEditor.
// /// @param database The path, including name, to the PHREEQC database file.
// auto setDatabase(std::string database) -> void;

// /// Set the aqueous species in the system by specifying which elements should exist.
// /// @param elements The names of the elements either as a vector of strings or as space-separated string list.
// auto setAqueousPhase(StringList elements) -> void;

// /// Set the gaseous phase in the system by specifying the end-member gases.
// /// @param elements The names of the gases either as a vector of strings or as space-separated string list.
// auto setGaseousPhase(StringList gases) -> void;

// /// Set the mineral phases in the system by specifying the names of the pure minerals.
// /// @param elements The names of the pure minerals either as a vector of strings or as space-separated string list.
// auto setMineralPhases(StringList minerals) -> void;

// /// Convert this PhreeqcEditor instance into a ChemicalSystem instance
// operator ChemicalSystem() const;

// /// Convert this PhreeqcEditor instance into a Phreeqc instance
// operator Phreeqc() const;

// private:
// struct Impl;

// std::unique_ptr<Impl> pimpl;
// };

// } // namespace Reaktoro
