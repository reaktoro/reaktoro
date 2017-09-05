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

/// The result of a thermodynamic model function that calculates standard thermodynamic properties of species.
template<typename VectorType>
struct ThermoModelResultBase
{
    /// Construct a default ThermoModelResultBase instance
    ThermoModelResultBase()
    {}

    /// Construct a ThermoModelResultBase instance with allocated memory
    explicit ThermoModelResultBase(Index nspecies)
    : standard_partial_molar_gibbs_energies(nspecies),
      standard_partial_molar_enthalpies(nspecies),
      standard_partial_molar_volumes(nspecies),
      standard_partial_molar_heat_capacities_cp(nspecies),
      standard_partial_molar_heat_capacities_cv(nspecies)
    {}

    /// Resize this ThermoModelResultBase with a given number of species
    auto resize(Index nspecies) -> void
    {
        standard_partial_molar_gibbs_energies.resize(nspecies);
        standard_partial_molar_enthalpies.resize(nspecies);
        standard_partial_molar_volumes.resize(nspecies);
        standard_partial_molar_heat_capacities_cp.resize(nspecies);
        standard_partial_molar_heat_capacities_cv.resize(nspecies);
    }

    /// Return a view of the thermodynamic properties of a phase.
    /// @param ispecies The index of the first species in the phase.
    /// @param nspecies The number of species in the phase.
    auto map(Index ispecies, Index nspecies) -> ThermoModelResultBase<VectorMap>
    {
        return {
            standard_partial_molar_gibbs_energies.map(ispecies, nspecies),
            standard_partial_molar_enthalpies.map(ispecies, nspecies),
            standard_partial_molar_volumes.map(ispecies, nspecies),
            standard_partial_molar_heat_capacities_cp.map(ispecies, nspecies),
            standard_partial_molar_heat_capacities_cv.map(ispecies, nspecies),
        };
    }

    /// Return a view of the thermodynamic properties of a phase.
    /// @param ispecies The index of the first species in the phase.
    /// @param nspecies The number of species in the phase.
    auto map(Index ispecies, Index nspecies) const -> ThermoModelResultBase<VectorConstMap>
    {
        return {
            standard_partial_molar_gibbs_energies.map(ispecies, nspecies),
            standard_partial_molar_enthalpies.map(ispecies, nspecies),
            standard_partial_molar_volumes.map(ispecies, nspecies),
            standard_partial_molar_heat_capacities_cp.map(ispecies, nspecies),
            standard_partial_molar_heat_capacities_cv.map(ispecies, nspecies),
        };
    }

    // Auxiliary type
    using ThermoVectorType = ThermoVectorBase<VectorType, VectorType, VectorType>;

    /// The standard partial molar Gibbs energies of the species (in units of J/mol).
    ThermoVectorType standard_partial_molar_gibbs_energies;

    /// The standard partial molar enthalpies of the species (in units of J/mol).
    ThermoVectorType standard_partial_molar_enthalpies;

    /// The standard partial molar volumes of the species (in units of m3/mol).
    ThermoVectorType standard_partial_molar_volumes;

    /// The standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    ThermoVectorType standard_partial_molar_heat_capacities_cp;

    /// The standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    ThermoVectorType standard_partial_molar_heat_capacities_cv;
};

/// The thermodynamic properties of the species in a chemical system.
using ThermoModelResult = ThermoModelResultBase<Vector>;

/// The thermodynamic properties of the species in a phase.
using PhaseThermoModelResult = ThermoModelResultBase<VectorMap>;

/// The thermodynamic properties of the species in a phase (constant).
using PhaseThermoModelResultConst = ThermoModelResultBase<VectorConstMap>;

/// The signature of the thermodynamic model function that calculates the standard thermodynamic properties of the species in a phase.
using PhaseThermoModel = std::function<void(PhaseThermoModelResult&, Temperature, Pressure)>;

/// The signature of the thermodynamic model function that calculates the standard thermodynamic properties of the species in a chemical system.
using ThermoModel = std::function<void(ThermoModelResult&, Temperature, Pressure)>;

} // namespace Reaktoro
