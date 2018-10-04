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

// Reaktoro includes
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Thermodynamics/Species/ThermoData.hpp>

namespace Reaktoro {

/// A type to describe the attributes of a mineral species
class MineralSpecies : public Species
{
public:
    /// Construct a default MineralSpecies instance
    MineralSpecies();

    /// Construct an MineralSpecies instance from a Species instance
    MineralSpecies(const Species& species);

    /// Set the thermodynamic data of the mineral species.
    auto setThermoData(const MineralSpeciesThermoData& thermo) -> void;

    /// Return the thermodynamic data of the mineral species.
    auto thermoData() const -> const MineralSpeciesThermoData&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
