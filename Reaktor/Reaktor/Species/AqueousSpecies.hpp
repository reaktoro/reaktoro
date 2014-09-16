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
#include <map>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Species/BaseSpecies.hpp>
#include <Reaktor/Species/ThermoParams.hpp>

namespace Reaktor {

/// A type for storing the parameters of the HKF equation of state for a aqueous species
struct AqueousThermoParamsHKF
{
    /// The apparent standard molal Gibbs free energy of formation of the species from its elements (in units of cal/mol)
    double Gf;

    /// The apparent standard molal enthalpy of formation of the species from its elements (in units of cal/mol)
    double Hf;

    /// The standard molal entropy of the species at reference temperature and pressure (in units of cal/(mol*K))
    double Sr;

    /// The coefficient a1 of the HKF equation of state of the aqueous species (in units of cal/(mol*bar))
    double a1;

    /// The coefficient a2 of the HKF equation of state of the aqueous species (in units of cal/mol)
    double a2;

    /// The coefficient a3 of the HKF equation of state of the aqueous species (in units of (cal*K)/(mol*bar))
    double a3;

    /// The coefficient a4 of the HKF equation of state of the aqueous species (in units of (cal*K)/mol)
    double a4;

    /// The coefficient c1 of the HKF equation of state of the aqueous species (in units of cal/(mol*K))
    double c1;

    /// The coefficient c2 of the HKF equation of state of the aqueous species (in units of (cal*K)/mol)
    double c2;

    /// The conventional Born coefficient of the aqueous species at reference temperature 298.15 K and pressure 1 bar (in units of cal/mol)
    double wref;
};

/// A type for storing the thermodynamic properties of an aqueous species
struct AqueousThermoParams
{
	/// The thermodynamic properties of an aqueous species as interpolated data
	ThermoParamsSpecies interpolated;

	/// The thermodynamic parameters of the HKF model for an aqueous species
	Optional<AqueousThermoParamsHKF> hkf;
};

/// A type to describe the attributes of an aqueous species
struct AqueousSpecies : public BaseSpecies
{
    /// The dissociation of a neutral aqueous species into charged species.
    /// For example, the dissociation of the aqueous species CaCl<sub>2</sub>(aq)
    /// produces 1 atom of Ca<sup>2+</sup> and 2 atoms of Cl<sup>-</sup>.
    std::map<std::string, double> dissociation;

    /// The thermodynamic parameters of the aqueous species
    AqueousThermoParams thermoparams;
};

} // namespace Reaktor
