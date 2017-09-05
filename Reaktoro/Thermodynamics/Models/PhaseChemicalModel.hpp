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
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>

namespace Reaktoro {

// Forward declarations
template<typename VectorType, typename MatrixType>
struct ChemicalModelResultBase;

/// The chemical properties of the species in a chemical system.
using ChemicalModelResult = ChemicalModelResultBase<Vector, Matrix>;

/// The chemical properties of the species in a phase.
using PhaseChemicalModelResult = ChemicalModelResultBase<VectorMap, MatrixMap>;

/// The chemical properties of the species in a phase (constant).
using PhaseChemicalModelResultConst = ChemicalModelResultBase<VectorConstMap, MatrixConstMap>;

/// The signature of the chemical model function that calculates the chemical properties of the species in a phase.
using PhaseChemicalModel = std::function<void(PhaseChemicalModelResult&, Temperature, Pressure, const Vector&)>;

/// The signature of the chemical model function that calculates the chemical properties of the species in a chemical system.
using ChemicalModel = std::function<void(ChemicalModelResult&, Temperature, Pressure, const Vector&)>;

/// The result of a chemical model function that calculates the chemical properties of species.
template<typename VectorType, typename MatrixType>
struct ChemicalModelResultBase
{
    /// Construct a default ChemicalModelResultBase instance
    ChemicalModelResultBase();

    /// Construct a ChemicalModelResultBase instance with allocated memory
    /// @param nphases The number of phases in the chemical system.
    /// @param nspecies The number of species in the chemical system.
    explicit ChemicalModelResultBase(Index nphases, Index nspecies)
    : ln_activity_coefficients(nspecies),
      ln_activity_constants(nspecies),
      ln_activities(nspecies),
      molar_volume(nspecies),
      residual_molar_gibbs_energy(nspecies),
      residual_molar_enthalpy(nspecies),
      residual_molar_heat_capacity_cp(nspecies),
      residual_molar_heat_capacity_cv(nspecies)
    {}

    /// Resize this ChemicalModelResultBase with a given number of species
    /// @param nphases The number of phases in the chemical system.
    /// @param nspecies The number of species in the chemical system.
    auto resize(Index nphases, Index nspecies) -> void
    {
        ln_activity_coefficients.resize(nspecies);
        ln_activity_constants.resize(nspecies);
        ln_activities.resize(nspecies);
        molar_volume.resize(nphases, nspecies);
        residual_molar_gibbs_energy.resize(nphases, nspecies);
        residual_molar_enthalpy.resize(nphases, nspecies);
        residual_molar_heat_capacity_cp.resize(nphases, nspecies);
        residual_molar_heat_capacity_cv.resize(nphases, nspecies);
    }

    /// Return a view of the chemical properties of a phase.
    /// @param iphase The index of the phase.
    /// @param ispecies The index of the first species in the phase.
    /// @param nspecies The number of species in the phase.
    auto map(Index iphase, Index ispecies, Index nspecies) -> PhaseChemicalModelResult
    {
        return {
            ln_activity_coefficients.map(ispecies, nspecies),
            ln_activity_constants.map(ispecies, nspecies),
            ln_activities.map(ispecies, nspecies),
            molar_volume.rowmap(iphase, ispecies, nspecies),
            residual_molar_gibbs_energy.rowmap(iphase, ispecies, nspecies),
            residual_molar_enthalpy.rowmap(iphase, ispecies, nspecies),
            residual_molar_heat_capacity_cp.rowmap(iphase, ispecies, nspecies),
            residual_molar_heat_capacity_cv.rowmap(iphase, ispecies, nspecies)
        };
    }

    /// Return a view of the chemical properties of a phase.
    /// @param iphase The index of the phase.
    /// @param ispecies The index of the first species in the phase.
    /// @param nspecies The number of species in the phase.
    auto map(Index iphase, Index ispecies, Index nspecies) const -> PhaseChemicalModelResultConst
    {
        return {
            ln_activity_coefficients.map(ispecies, nspecies),
            ln_activity_constants.map(ispecies, nspecies),
            ln_activities.map(ispecies, nspecies),
            molar_volume.rowmap(iphase, ispecies, nspecies),
            residual_molar_gibbs_energy.rowmap(iphase, ispecies, nspecies),
            residual_molar_enthalpy.rowmap(iphase, ispecies, nspecies),
            residual_molar_heat_capacity_cp.rowmap(iphase, ispecies, nspecies),
            residual_molar_heat_capacity_cv.rowmap(iphase, ispecies, nspecies)
        };
    }

    // Auxiliary types
    using ThermoVectorType = ThermoVectorBase<VectorType, VectorType, VectorType>;
    using ChemicalVectorType = ChemicalVectorBase<VectorType, VectorType, VectorType, MatrixType>;

    /// The natural log of the activity coefficients of the species.
    ChemicalVectorType ln_activity_coefficients;

    /// The natural log of the activity constants of the species.
    ThermoVectorType ln_activity_constants;

    /// The natural log of the activities of the species.
    ChemicalVectorType ln_activities;

    /// The molar volumes of the phases (in units of m3/mol).
    ChemicalVectorType molar_volume; // TODO change molar_volume to phase_molar_volumes

    /// The residual molar Gibbs energies of the phases w.r.t. to its ideal state (in units of J/mol).
    ChemicalVectorType residual_molar_gibbs_energy; // TODO change residual_molar_gibbs_energy to phase_residual_molar_gibbs_energies

    /// The residual molar enthalpies of the phases w.r.t. to its ideal state (in units of J/mol).
    ChemicalVectorType residual_molar_enthalpy; // TODO change residual_molar_enthalpy to phase_residual_molar_enthalpies

    /// The residual molar isobaric heat capacities of the phases w.r.t. to its ideal state (in units of J/(mol*K)).
    ChemicalVectorType residual_molar_heat_capacity_cp; // TODO change residual_molar_heat_capacity_cp to phase_residual_molar_heat_capacities_cp

    /// The residual molar isochoric heat capacities of the phases w.r.t. to its ideal state (in units of J/(mol*K)).
    ChemicalVectorType residual_molar_heat_capacity_cv; // TODO change residual_molar_heat_capacity_cv to phase_residual_molar_heat_capacities_cv
};

} // namespace Reaktoro
