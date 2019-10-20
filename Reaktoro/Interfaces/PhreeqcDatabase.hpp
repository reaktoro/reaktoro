// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <memory>
#include <set>
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

// Forward declarations
class Database;
class Element;
class AqueousSpecies;
class FluidSpecies;
using GaseousSpecies = FluidSpecies;
class MineralSpecies;

class PhreeqcDatabase
{
public:
    /// Construct a default PhreeqcDatabase instance
    PhreeqcDatabase();

    /// Construct a custom PhreeqcDatabase instance
    /// @param filename The path to the Phreeqc database file
    explicit PhreeqcDatabase(std::string filename);

    /// Load a Phreeqc database.
    /// @param filename The path to the Phreeqc database file
    auto load(std::string filename) -> void;

    auto numElements() const -> unsigned;

    auto numAqueousSpecies() const -> unsigned;

    auto numGaseousSpecies() const -> unsigned;

    auto numMineralSpecies() const -> unsigned;

    auto numMasterSpecies() const -> unsigned;

    auto numProductSpecies() const -> unsigned;

    auto element(Index index) const -> Element;

    auto elements() const -> const std::vector<Element>&;

    auto aqueousSpecies(Index index) const -> AqueousSpecies;

    auto aqueousSpecies(std::string name) const -> AqueousSpecies;

    auto aqueousSpecies() const -> const std::vector<AqueousSpecies>&;

    auto gaseousSpecies(Index index) const -> GaseousSpecies;

    auto gaseousSpecies(std::string name) const -> GaseousSpecies;

    auto gaseousSpecies() const -> const std::vector<GaseousSpecies>&;

    auto mineralSpecies(Index index) const -> MineralSpecies;

    auto mineralSpecies(std::string name) const -> MineralSpecies;

    auto containsAqueousSpecies(std::string name) const -> bool;

    auto containsGaseousSpecies(std::string name) const -> bool;

    auto containsMineralSpecies(std::string name) const -> bool;

    auto mineralSpecies() const -> const std::vector<MineralSpecies>&;

    auto masterSpecies() const -> std::set<std::string>;

    /// Cross this PhreeqcDatabase instance with master thermodynamic data in another Database instance
    auto cross(const Database& master) -> Database;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
