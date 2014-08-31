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

#include "ThermoData.hpp"

namespace Reaktor {

// Forward declarations
class Database;
class AqueousSpecies;
class GaseousSpecies;
class MineralSpecies;
class BilinearInterpolator;
struct ThermoDataSpecies;

class ThermoEditor
{
public:
    explicit ThermoEditor(const Database& database);

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

    ThermoDataSpecies
    thermoDataHKF(const AqueousSpecies& species) const;

    ThermoDataSpecies
    thermoDataHKF(const GaseousSpecies& species) const;

    ThermoDataSpecies
    thermoDataHKF(const MineralSpecies& species) const;

    ThermoDataSpecies
    thermoDataReaction(const AqueousSpecies& species) const;

    ThermoDataSpecies
    thermoDataReaction(const GaseousSpecies& species) const;

    ThermoDataSpecies
    thermoDataReaction(const MineralSpecies& species) const;

    ThermoDataSpecies
    thermoDataInterpolated(const AqueousSpecies& species) const;

    ThermoDataSpecies
    thermoDataInterpolated(const GaseousSpecies& species) const;

    ThermoDataSpecies
    thermoDataInterpolated(const MineralSpecies& species) const;

    ThermoDataSpecies
    thermoData(const std::string& species) const;

    ThermoDataSpecies
    thermoData(const AqueousSpecies& species) const;

    ThermoDataSpecies
    thermoData(const GaseousSpecies& species) const;

    ThermoDataSpecies
    thermoData(const MineralSpecies& species) const;

private:
    /// The database instance
    const Database& database$;

    /// The temperatures for the interpolation (in units of K)
    std::vector<double> temperatures$;

    /// The pressures for the interpolation (in units of Pa)
    std::vector<double> pressures$;
};


} // namespace Reaktor
