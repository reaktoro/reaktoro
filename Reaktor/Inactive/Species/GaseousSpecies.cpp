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

#include "GaseousSpecies.hpp"

// Reaktor includes
#include <Reaktor/Math/BilinearInterpolator.hpp>
#include <Reaktor/Thermo/ThermoUtils.hpp>

namespace Reaktor {

GaseousSpecies::GaseousSpecies()
{}

auto GaseousSpecies::setGas(const std::string& gas) -> void
{
    gas$ = gas;
}

auto GaseousSpecies::setThermoData(const GaseousThermoData& thermoData) -> void
{
    thermo_data$ = thermoData;
}

auto GaseousSpecies::gas() const -> const std::string&
{
    return gas$;
}

auto GaseousSpecies::thermoData() const -> const GaseousThermoData&
{
    return thermo_data$;
}

auto operator<<(std::ostream& out, const GaseousSpecies& species) -> std::ostream&
{
    // Get the HKF thermodynamic data of the species
    const GaseousThermoDataHKF& hkf = species.thermoData().hkf.get();

    out << static_cast<GeneralSpecies>(species);
    out << "  gas name: " << species.gas() << std::endl;
    out << "  thermo data (HKF)" << std::endl;
    out << "    Gf: " << hkf.Gf << std::endl;
    out << "    Hf: " << hkf.Hf << std::endl;
    out << "    Sr: " << hkf.Sr << std::endl;
    out << "    a: " << hkf.a << std::endl;
    out << "    b: " << hkf.b << std::endl;
    out << "    c: " << hkf.c << std::endl;
    out << "    Tmax: " << hkf.Tmax << std::endl;

    return out;
}

} // namespace Reaktor

