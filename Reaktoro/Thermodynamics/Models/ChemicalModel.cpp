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

#include "ChemicalModel.hpp"

namespace Reaktoro {

ChemicalModelResult::ChemicalModelResult()
{}

ChemicalModelResult::ChemicalModelResult(Index nphases, Index nspecies)
: ln_activity_coefficients(nspecies),
  ln_activities(nspecies),
  phase_molar_volumes(nspecies),
  phase_residual_molar_gibbs_energies(nspecies),
  phase_residual_molar_enthalpies(nspecies),
  phase_residual_molar_heat_capacities_cp(nspecies),
  phase_residual_molar_heat_capacities_cv(nspecies)
{}

auto ChemicalModelResult::resize(Index nphases, Index nspecies) -> void
{
    ln_activity_coefficients.resize(nspecies);
    ln_activities.resize(nspecies);
    phase_molar_volumes.resize(nphases, nspecies);
    phase_residual_molar_gibbs_energies.resize(nphases, nspecies);
    phase_residual_molar_enthalpies.resize(nphases, nspecies);
    phase_residual_molar_heat_capacities_cp.resize(nphases, nspecies);
    phase_residual_molar_heat_capacities_cv.resize(nphases, nspecies);
}

auto ChemicalModelResult::phaseProperties(Index iphase, Index ispecies, Index nspecies) -> PhaseChemicalModelResult
{
    return {
        rows(ln_activity_coefficients, ispecies, nspecies),
        rows(ln_activities, ispecies, nspecies),
        row(phase_molar_volumes, iphase, ispecies, nspecies),
        row(phase_residual_molar_gibbs_energies, iphase, ispecies, nspecies),
        row(phase_residual_molar_enthalpies, iphase, ispecies, nspecies),
        row(phase_residual_molar_heat_capacities_cp, iphase, ispecies, nspecies),
        row(phase_residual_molar_heat_capacities_cv, iphase, ispecies, nspecies)
    };
}

auto ChemicalModelResult::phaseProperties(Index iphase, Index ispecies, Index nspecies) const -> PhaseChemicalModelResultConst
{
    return {
        rows(ln_activity_coefficients, ispecies, nspecies),
        rows(ln_activities, ispecies, nspecies),
        row(phase_molar_volumes, iphase, ispecies, nspecies),
        row(phase_residual_molar_gibbs_energies, iphase, ispecies, nspecies),
        row(phase_residual_molar_enthalpies, iphase, ispecies, nspecies),
        row(phase_residual_molar_heat_capacities_cp, iphase, ispecies, nspecies),
        row(phase_residual_molar_heat_capacities_cv, iphase, ispecies, nspecies)
    };
}

} // namespace Reaktoro
