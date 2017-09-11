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

#include "ThermoModel.hpp"

namespace Reaktoro {

ThermoModelResult::ThermoModelResult()
{}

ThermoModelResult::ThermoModelResult(Index nspecies)
: standard_partial_molar_gibbs_energies(nspecies),
  standard_partial_molar_enthalpies(nspecies),
  standard_partial_molar_volumes(nspecies),
  standard_partial_molar_heat_capacities_cp(nspecies),
  standard_partial_molar_heat_capacities_cv(nspecies),
  ln_activity_constants(nspecies)
{}

auto ThermoModelResult::resize(Index nspecies) -> void
{
    standard_partial_molar_gibbs_energies.resize(nspecies);
    standard_partial_molar_enthalpies.resize(nspecies);
    standard_partial_molar_volumes.resize(nspecies);
    standard_partial_molar_heat_capacities_cp.resize(nspecies);
    standard_partial_molar_heat_capacities_cv.resize(nspecies);
    ln_activity_constants.resize(nspecies);
}

auto ThermoModelResult::phaseProperties(Index iphase, Index ispecies, Index nspecies) -> PhaseThermoModelResult
{
    return {
        rows(standard_partial_molar_gibbs_energies, ispecies, nspecies),
        rows(standard_partial_molar_enthalpies, ispecies, nspecies),
        rows(standard_partial_molar_volumes, ispecies, nspecies),
        rows(standard_partial_molar_heat_capacities_cp, ispecies, nspecies),
        rows(standard_partial_molar_heat_capacities_cv, ispecies, nspecies),
        rows(ln_activity_constants, ispecies, nspecies),
    };
}

auto ThermoModelResult::phaseProperties(Index iphase, Index ispecies, Index nspecies) const -> PhaseThermoModelResultConst
{
    return {
        rows(standard_partial_molar_gibbs_energies, ispecies, nspecies),
        rows(standard_partial_molar_enthalpies, ispecies, nspecies),
        rows(standard_partial_molar_volumes, ispecies, nspecies),
        rows(standard_partial_molar_heat_capacities_cp, ispecies, nspecies),
        rows(standard_partial_molar_heat_capacities_cv, ispecies, nspecies),
        rows(ln_activity_constants, ispecies, nspecies),
    };
}

} // namespace Reaktoro
