// Reaktoro is a C++ library for computational reaction modelling.
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
#include <functional>

// Reaktoro includes
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

/// The result of the thermodynamic model function that calculates the standard thermodynamic properties of a phase.
struct PhaseThermoModelResult
{
    /// Construct a default PhaseThermoModelResult instance
    PhaseThermoModelResult();

    /// Construct a PhaseThermoModelResult instance with allocated memory
    explicit PhaseThermoModelResult(unsigned nspecies);

    /// Resize this PhaseThermoModelResult with a given number of species
    auto resize(unsigned nspecies) -> void;

    /// The standard partial molar Gibbs energies of the species (in units of J/mol).
    ThermoVector standard_partial_molar_gibbs_energies;

    /// The standard partial molar enthalpies of the species (in units of J/mol).
    ThermoVector standard_partial_molar_enthalpies;

    /// The standard partial molar volumes of the species (in units of m3/mol).
    ThermoVector standard_partial_molar_volumes;

    /// The standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    ThermoVector standard_partial_molar_heat_capacities_cp;

    /// The standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    ThermoVector standard_partial_molar_heat_capacities_cv;
};

/// The signature of the thermodynamic model function that calculates the standard thermodynamic properties of a phase.
using PhaseThermoModel = std::function<PhaseThermoModelResult(double, double)>;

} // namespace Reaktoro
