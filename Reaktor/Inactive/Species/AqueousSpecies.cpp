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

#include "AqueousSpecies.hpp"

// C++ includes
#include <sstream>

// Reaktor includes
#include <Reaktor/Math/BilinearInterpolator.hpp>
#include <Reaktor/Thermo/ThermoUtils.hpp>
#include <Reaktor/Utils/StringUtils.hpp>

namespace Reaktor {

AqueousSpecies::AqueousSpecies()
{}

auto AqueousSpecies::setDissociation(const std::map<std::string, double>& dissociation) -> void
{
    dissociation$ = dissociation;
}

auto AqueousSpecies::setDissociation(const std::string& dissociation) -> void
{
    // Reset the dissociation
    dissociation$.clear();

    // Split the complex formula in words delimited by ',' and ' '
    auto words = split(dissociation, ", ");

    // Create the pair entries
    for(const auto& word : words)
    {
        auto pair = split(word, ":");
        dissociation$.insert({pair[1], tofloat(pair[0])});
    }
}

auto AqueousSpecies::setThermoData(const AqueousThermoData& thermoData) -> void
{
    thermo_data$ = thermoData;
}

auto AqueousSpecies::thermoData() const -> const AqueousThermoData&
{
    return thermo_data$;
}

auto AqueousSpecies::dissociation() const -> const std::map<std::string, double>&
{
    return dissociation$;
}

auto operator<<(std::ostream& out, const AqueousSpecies& species) -> std::ostream&
{
    auto dissociationstr = [](const AqueousSpecies& species)
    {
        std::stringstream ss;
        for(const auto& pair : species.dissociation())
            ss << pair.second << ":" << pair.second << ", ";
        return ss.str().substr(0, ss.str().size() - 2);
    };

    // Get the HKF thermodynamic data of the species
    const AqueousThermoDataHKF& hkf = species.thermoData().hkf.get();

    out << static_cast<GeneralSpecies>(species);
    out << "  electrical charge: " << species.charge() << std::endl;
    out << "  dissociation: " << dissociationstr(species) << std::endl;
    out << "  thermo data (HKF)" << std::endl;
    out << "    Gf: " << hkf.Gf << std::endl;
    out << "    Hf: " << hkf.Hf << std::endl;
    out << "    Sr: " << hkf.Sr << std::endl;
    out << "    a1: " << hkf.a1 << std::endl;
    out << "    a2: " << hkf.a2 << std::endl;
    out << "    a3: " << hkf.a3 << std::endl;
    out << "    a4: " << hkf.a4 << std::endl;
    out << "    c1: " << hkf.c1 << std::endl;
    out << "    c2: " << hkf.c2 << std::endl;
    out << "    wref: " << hkf.wref << std::endl;
    return out;
}

} // namespace Reaktor
