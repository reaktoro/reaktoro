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

// Reaktoro includes
#include "GeneralSpecies.hpp"
#include "ThermoData.hpp"

namespace Reaktoro {

/// A type to describe the attributes of a mineral species
class MineralSpecies : public GeneralSpecies
{
public:
    /// Construct a default MineralSpecies instance
    MineralSpecies();

    /// Construct an MineralSpecies instance from a Species instance
    MineralSpecies(const GeneralSpecies& species);

    /// Construct a copy of an MineralSpecies instance
    MineralSpecies(const MineralSpecies& other);

    /// Destroy this instance
    virtual ~MineralSpecies();

    /// Assign an MineralSpecies instance to this instance
    auto operator=(MineralSpecies other) -> MineralSpecies&;

    /// Set the thermodynamic data of the mineral species.
    auto setThermoData(const MineralSpeciesThermoData& thermo) -> void;

    /// Return the thermodynamic data of the mineral species.
    auto thermoData() const -> const MineralSpeciesThermoData&;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
