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
#include <vector>
#include <string>

// Reaktor includes
#include <Reaktor/Thermo/ThermoEditor.hpp>

namespace Reaktor {

// Reaktor forward declarations
class Phase;
class Species;
class Database;
class AqueousPhase;
class GaseousPhase;
class MineralPhase;
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;
class BilinearInterpolator;

class ChemicalConverter
{
public:
    explicit ChemicalConverter(const Database& database);

    void
    setTemperatures(const std::vector<double>& temperatures);

    void
    setPressures(const std::vector<double>& pressures);

    const Database&
    database() const;

    const std::vector<double>&
    temperatures() const;

    const std::vector<double>&
    pressures() const;

    Species
    convertSpecies(const std::string& name) const;

    /**
     * Converts an AqueousSpecies instance to a Species instance
     * @param species The aqueous species instance
     * @return A Species instance created from the given aqueous species
     */
    Species
    convertSpecies(const AqueousSpecies& species) const;

    /**
     * Converts a GaseousSpecies instance to a Species instance
     * @param species The gaseous species instance
     * @return A Species instance created from the given gaseous species
     */
    Species
    convertSpecies(const GaseousSpecies& species) const;

    /**
     * Converts a MineralSpecies instance to a Species instance
     * @param species The mineral species instance
     * @return A Species instance created from the given mineral species
     */
    Species
    convertSpecies(const MineralSpecies& species) const;

    /**
     * Converts an AqueousPhase instance to a Phase instance
     * @param phase The aqueous phase instance
     * @return A Phase instance created from the given aqueous phase
     */
    Phase
    convertPhase(const AqueousPhase& phase) const;

    /**
     * Converts a GaseousPhase instance to a Phase instance
     * @param phase The gaseous phase instance
     * @return A Phase instance created from the given gaseous phase
     */
    Phase
    convertPhase(const GaseousPhase& phase) const;

    /**
     * Converts a MineralPhase instance to a Phase instance
     * @param phase The mineral phase instance
     * @return A Phase instance created from the given mineral phase
     */
    Phase
    convertPhase(const MineralPhase& phase) const;

private:
    ThermoEditor thermo_editor$;
};

} /* namespace Reaktor */
