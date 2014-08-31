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

#include "ThermoStateSpecies.hpp"

// Reaktor includes
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>
#include <Reaktor/Thermodynamics/ThermoStateSpeciesHKF.hpp>
#include <Reaktor/Thermodynamics/WaterConstants.hpp>
#include <Reaktor/Thermodynamics/ThermoStateWater.hpp>

namespace Reaktor
{
ThermoStateSpecies::ThermoStateSpecies()
: volume(0), entropy(0), helmholtz(0), internal_energy(0), enthalpy(0), gibbs(0), cp(0)
{}

auto operator<<(std::ostream& out, const ThermoStateSpecies& st) -> std::ostream&
{
	out << "volume          = " << st.volume          << std::endl;
	out << "entropy         = " << st.entropy         << std::endl;
	out << "helmholtz       = " << st.helmholtz       << std::endl;
	out << "internal_energy = " << st.internal_energy << std::endl;
	out << "enthalpy        = " << st.enthalpy        << std::endl;
	out << "gibbs           = " << st.gibbs           << std::endl;
	out << "cp              = " << st.cp              << std::endl;

	return out;
}

auto speciesThermo(double T, double P, const AqueousSpecies& species) -> ThermoStateSpecies
{
    // Check if the given (T, P) falls inside the gaseous region
    if(T < waterCriticalTemperature)
    {
        // Calculate the saturated vapor pressure
        const double Psat = saturatedPressureWater(T);

        // Set the pressure to one pascal higher than the saturated vapor pressure
        P = (P < Psat) ? Psat + 1 : P;
    }

    return speciesThermoHKF(T, P, species);
}

auto speciesThermo(double T, double P, const GaseousSpecies& species) -> ThermoStateSpecies
{
    return speciesThermoHKF(T, P, species);
}

auto speciesThermo(double T, double P, const MineralSpecies& species) -> ThermoStateSpecies
{
    return speciesThermoHKF(T, P, species);
}

} // namespace Reaktor
