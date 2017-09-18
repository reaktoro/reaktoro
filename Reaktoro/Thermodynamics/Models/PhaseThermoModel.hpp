// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

/// The result of a thermodynamic model function that calculates the thermodynamic properties of species.
template<typename VectorType>
struct PhaseThermoModelResultBase
{
    /// The standard partial molar Gibbs energies of the species (in units of J/mol).
    VectorType standard_partial_molar_gibbs_energies;

    /// The standard partial molar enthalpies of the species (in units of J/mol).
    VectorType standard_partial_molar_enthalpies;

    /// The standard partial molar volumes of the species (in units of m3/mol).
    VectorType standard_partial_molar_volumes;

    /// The standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    VectorType standard_partial_molar_heat_capacities_cp;

    /// The standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    VectorType standard_partial_molar_heat_capacities_cv;

    /// The natural log of the activity constants of the species.
    VectorType ln_activity_constants;
};

/// The thermodynamic properties of the species in a phase.
using PhaseThermoModelResult = PhaseThermoModelResultBase<ThermoVectorRef>;

/// The thermodynamic properties of the species in a phase (constant).
using PhaseThermoModelResultConst = PhaseThermoModelResultBase<ThermoVectorConstRef>;

/// The signature of the chemical model function that calculates the thermodynamic properties of the species in a phase.
using PhaseThermoModel = std::function<void(PhaseThermoModelResult&, Temperature, Pressure)>;

} // namespace Reaktoro
