// Reaktor is a C++ library for computational reaction modelling.
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

// Reaktor includes
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Species/GeneralSpecies.hpp>
#include <Reaktor/Thermodynamics/ThermoData.hpp>

namespace Reaktor {

class MineralSpecies : public GeneralSpecies
{
public:
    MineralSpecies();

    auto setMolarVolume(units::MolarVolume value) -> void;

    auto setThermoData(const MineralThermoData& thermoData) -> void;

    auto density() const -> units::Density;

    auto molarVolume() const -> units::MolarVolume;

    auto thermoData() const -> const MineralThermoData&;

private:
    /// The molar volume of the mineral
    units::MolarVolume molar_volume$;

    /// The mineral thermodynamic data of the species from the HKF model
    MineralThermoData thermo_data$;
};

/**
 * Outputs a MineralSpecies instance
 */
auto operator<<(std::ostream& out, const MineralSpecies& species) -> std::ostream&;

} // namespace Reaktor
