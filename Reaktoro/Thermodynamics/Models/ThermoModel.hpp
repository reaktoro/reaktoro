// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <functional>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/PhaseThermoModel.hpp>

namespace Reaktoro {

/// The result of a thermodynamic model function that calculates standard thermodynamic properties of species.
class ThermoModelResult
{
public:
    /// Construct a default ThermoModelResult instance.
    ThermoModelResult();

    /// Construct a ThermoModelResult instance with allocated memory.
    /// @param nphases The number of phases in the chemical system.
    /// @param nspecies The number of species in the chemical system.
    ThermoModelResult(Index nphases, Index nspecies);

    /// Resize this ThermoModelResult with a given number of species
    /// @param nphases The number of phases in the chemical system.
    /// @param nspecies The number of species in the chemical system.
    auto resize(Index nphases, Index nspecies) -> void;

    /// Return a view of the thermodynamic properties of a phase.
    /// @param iphase The index of the phase.
    /// @param ispecies The index of the first species in the phase.
    /// @param nspecies The number of species in the phase.
    auto phaseProperties(Index iphase, Index ispecies, Index nspecies) -> PhaseThermoModelResult;

    /// Return a view of the thermodynamic properties of a phase.
    /// @param iphase The index of the phase.
    /// @param ispecies The index of the first species in the phase.
    /// @param nspecies The number of species in the phase.
    auto phaseProperties(Index iphase, Index ispecies, Index nspecies) const -> PhaseThermoModelResultConst;

    /// The standard partial molar Gibbs energies of the species (in units of J/mol).
    inline auto standardPartialMolarGibbsEnergies() -> ThermoVectorRef { return standard_partial_molar_gibbs_energies; };

    /// The standard partial molar Gibbs energies of the species (in units of J/mol).
    inline auto standardPartialMolarGibbsEnergies() const -> ThermoVectorConstRef { return standard_partial_molar_gibbs_energies; };

    /// The standard partial molar enthalpies of the species (in units of J/mol).
    inline auto standardPartialMolarEnthalpies() -> ThermoVectorRef { return standard_partial_molar_enthalpies; };

    /// The standard partial molar enthalpies of the species (in units of J/mol).
    inline auto standardPartialMolarEnthalpies() const -> ThermoVectorConstRef { return standard_partial_molar_enthalpies; };

    /// The standard partial molar volumes of the species (in units of m3/mol).
    inline auto standardPartialMolarVolumes() -> ThermoVectorRef { return standard_partial_molar_volumes; };

    /// The standard partial molar volumes of the species (in units of m3/mol).
    inline auto standardPartialMolarVolumes() const -> ThermoVectorConstRef { return standard_partial_molar_volumes; };

    /// The standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    inline auto standardPartialMolarHeatCapacitiesConstP() -> ThermoVectorRef { return standard_partial_molar_heat_capacities_cp; };

    /// The standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    inline auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVectorConstRef { return standard_partial_molar_heat_capacities_cp; };

    /// The standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    inline auto standardPartialMolarHeatCapacitiesConstV() -> ThermoVectorRef { return standard_partial_molar_heat_capacities_cv; };

    /// The standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    inline auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVectorConstRef { return standard_partial_molar_heat_capacities_cv; };

    /// The natural log of the activity constants of the species.
    inline auto lnActivityConstants() -> ThermoVectorRef { return ln_activity_constants; };

    /// The natural log of the activity constants of the species.
    inline auto lnActivityConstants() const -> ThermoVectorConstRef { return ln_activity_constants; };

private:
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

    /// The natural log of the activity constants of the species.
    ThermoVector ln_activity_constants;
};

/// The signature of the thermodynamic model function that calculates the standard thermodynamic properties of the species in a chemical system.
using ThermoModel = std::function<void(ThermoModelResult&, Temperature, Pressure)>;

} // namespace Reaktoro
