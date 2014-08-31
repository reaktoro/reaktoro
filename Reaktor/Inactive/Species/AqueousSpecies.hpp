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

// C++ includes
#include <map>
#include <string>

// Reaktor includes
#include <Reaktor/Species/GeneralSpecies.hpp>
#include <Reaktor/Thermo/ThermoData.hpp>

namespace Reaktor {

class AqueousSpecies : public GeneralSpecies
{
public:
    AqueousSpecies();

    /**
     * Sets the ionic dissociation of the aqueous species
     *
     * The ionic dissociation of an aqueous species is described as a
     * sequence of pairs (@e ion, @e stoichiometry). It is shown below
     * how the ionic dissociation of species NaCl(aq) can be set using
     * this method:
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * AqueousSpecies species;
     * species.setName("NaCl(aq)");
     * species.setDissociation({{"Na+", 1}, {"Cl-", 1}});
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param dissociation The ionic dissociation of the aqueous species
     */
    auto setDissociation(const std::map<std::string, double>& dissociation) -> void;

    /**
     * Sets the ionic dissociation of the aqueous species
     *
     * The ionic dissociation of an aqueous species is described as a
     * sequence of pairs (@e ion, @e stoichiometry). It is shown below
     * how the ionic dissociation of species NaCl(aq) can be set using
     * this method:
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * AqueousSpecies species;
     * species.setName("NaCl(aq)");
     * species.setDissociation("1:Na+, 1:Cl-");
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param dissociation The ionic dissociation of the aqueous species as a formated string
     */
    auto setDissociation(const std::string& dissociation) -> void;

    auto setThermoData(const AqueousThermoData& thermoData) -> void;

    auto thermoData() const -> const AqueousThermoData&;

    auto dissociation() const -> const std::map<std::string, double>&;

private:
    /// The dissociation formula of the aqueous species in case it is an aqueous complex
    std::map<std::string, double> dissociation$;

    /// The thermodynamic data of the aqueous species from the HKF model
    AqueousThermoData thermo_data$;
};

/**
 * Outputs an AqueousSpecies instance
 */
auto operator<<(std::ostream& out, const AqueousSpecies& species) -> std::ostream&;

} // namespace Reaktor
