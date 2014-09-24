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

// C++ includes
#include <string>

// Reaktor includes
#include <Reaktor/Species/BaseSpecies.hpp>
#include <Reaktor/Species/ThermoParams.hpp>

namespace Reaktor {

/// A type for storing the parameters of the HKF equation of state for a gaseous species
struct GaseousThermoParamsHKF
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    double Gf;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    double Hf;

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol�K))
    double Sr;

    /// The coefficient a of the HKF equation of state of the gaseous species (in units of cal/(mol�K))
    double a;

    /// The coefficient b of the HKF equation of state of the gaseous species (in units of cal/(mol�K^2))
    double b;

    /// The coefficient c of the HKF equation of state of the gaseous species (in units of (cal�K)/mol)
    double c;

    /// The maximum temperature at which the HKF equation of state can be applied for the gaseous species (in units of K)
    double Tmax;
};

/// A type for storing the thermodynamic properties of a gaseous species
struct GaseousThermoParams
{
	/// The thermodynamic properties of a gaseous species as interpolated data
	ThermoParamsSpecies interpolated;

	/// The thermodynamic parameters of the HKF model for a gaseous species
	Optional<GaseousThermoParamsHKF> hkf;
};

/// A type to describe the attributes of a gaseous species
struct GaseousSpecies : public BaseSpecies
{
    /// The technical name of the gas
    std::string gas;

    /// The thermodynamic parameters of the gaseous species
    GaseousThermoParams thermoparams;
};

} // namespace Reaktor
