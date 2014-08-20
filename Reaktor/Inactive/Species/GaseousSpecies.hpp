/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// Reaktor includes
#include <Reaktor/Species/GeneralSpecies.hpp>
#include <Reaktor/Thermo/ThermoData.hpp>

namespace Reaktor {

class GaseousSpecies : public GeneralSpecies
{
public:
    GaseousSpecies();

    auto setGas(const std::string& gas) -> void;

    auto setThermoData(const GaseousThermoData& thermoData) -> void;

    auto gas() const -> const std::string&;

    auto thermoData() const -> const GaseousThermoData&;

private:
    /// The technical name of the gas
    std::string gas$;

    /// The thermodynamic data of the gaseous species
    GaseousThermoData thermo_data$;
};

/**
 * Outputs a GaseousSpecies instance
 */
auto operator<<(std::ostream& out, const GaseousSpecies& species) -> std::ostream&;

} /* namespace Reaktor */
