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

#include "PhaseChemicalModel.hpp"

namespace Reaktoro {

PhaseChemicalModelResult::PhaseChemicalModelResult()
{}

PhaseChemicalModelResult::PhaseChemicalModelResult(unsigned nspecies)
: ln_activity_coefficients(nspecies),
  ln_activities(nspecies),
  molar_volume(nspecies),
  residual_molar_gibbs_energy(nspecies),
  residual_molar_enthalpy(nspecies),
  residual_molar_heat_capacity_cp(nspecies),
  residual_molar_heat_capacity_cv(nspecies)
{}

auto PhaseChemicalModelResult::resize(unsigned nspecies) -> void
{
    ln_activity_coefficients.resize(nspecies);
    ln_activities.resize(nspecies);
    molar_volume = ChemicalScalar(nspecies);
    residual_molar_gibbs_energy = ChemicalScalar(nspecies);
    residual_molar_enthalpy = ChemicalScalar(nspecies);
    residual_molar_heat_capacity_cp = ChemicalScalar(nspecies);
    residual_molar_heat_capacity_cv = ChemicalScalar(nspecies);
}

} // namespace Reaktoro
