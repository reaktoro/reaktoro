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
#include <map>
#include <string>

// Reaktoro includes
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Thermodynamics/Species/ThermoData.hpp>

namespace Reaktoro {

/// A type to represent an aqueous species
class AqueousSpecies : public Species
{
public:
    /// Construct a default AqueousSpecies instance
    AqueousSpecies();

    /// Construct an AqueousSpecies instance from a Species instance
    AqueousSpecies(const Species& species);

    /// Set the charge of the aqueous species.
    auto setCharge(double value) -> void;

    /// Set the dissociation of a neutral aqueous species into charged species.
    /// For example, the dissociation of the aqueous species CaCl<sub>2</sub>(aq)
    /// produces 1 atom of Ca<sup>2+</sup> and 2 atoms of Cl<sup>-</sup>.
    auto setDissociation(const std::map<std::string, double>& dissociation) -> void;

    /// Set the thermodynamic data of the aqueous species.
    auto setThermoData(const AqueousSpeciesThermoData& thermo) -> void;

    /// Return the electrical charge of the aqueous species
    auto charge() const -> double;

    /// Return the dissociation of the aqueous species.
    auto dissociation() const -> const std::map<std::string, double>&;

    /// Return the thermodynamic data of the aqueous species.
    auto thermoData() const -> const AqueousSpeciesThermoData&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
