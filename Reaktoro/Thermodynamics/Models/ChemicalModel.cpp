// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2020 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #include "ChemicalModel.hpp"

// namespace Reaktoro {

// ChemicalModelResult::ChemicalModelResult()
// {}

// ChemicalModelResult::ChemicalModelResult(Index nphases, Index nspecies)
// : ln_activity_coefficients(nspecies),
//   ln_activities(nspecies),
//   partial_molar_volumes(nspecies),
//   phase_molar_volumes(nphases, nspecies),
//   phase_residual_molar_gibbs_energies(nphases, nspecies),
//   phase_residual_molar_enthalpies(nphases, nspecies),
//   phase_residual_molar_heat_capacities_cp(nphases, nspecies),
//   phase_residual_molar_heat_capacities_cv(nphases, nspecies)
// {}

// auto ChemicalModelResult::resize(Index nphases, Index nspecies) -> void
// {
//     ln_activity_coefficients.resize(nspecies);
//     ln_activities.resize(nspecies);
//     partial_molar_volumes.resize(nspecies);
//     phase_molar_volumes.resize(nphases, nspecies);
//     phase_residual_molar_gibbs_energies.resize(nphases, nspecies);
//     phase_residual_molar_enthalpies.resize(nphases, nspecies);
//     phase_residual_molar_heat_capacities_cp.resize(nphases, nspecies);
//     phase_residual_molar_heat_capacities_cv.resize(nphases, nspecies);
// }

// auto ChemicalModelResult::phaseProperties(Index iphase, Index ispecies, Index nspecies)-> ActivityModelFnResult
// {
//     return {
//         rows(ln_activity_coefficients, ispecies, ispecies, nspecies, nspecies),
//         rows(ln_activities, ispecies, ispecies, nspecies, nspecies),
//         rows(partial_molar_volumes, ispecies, ispecies, nspecies, nspecies),
//         row(phase_molar_volumes, iphase, ispecies, nspecies),
//         row(phase_residual_molar_gibbs_energies, iphase, ispecies, nspecies),
//         row(phase_residual_molar_enthalpies, iphase, ispecies, nspecies),
//         row(phase_residual_molar_heat_capacities_cp, iphase, ispecies, nspecies),
//         row(phase_residual_molar_heat_capacities_cv, iphase, ispecies, nspecies)
//     };
// }

// auto ChemicalModelResult::phaseProperties(Index iphase, Index ispecies, Index nspecies) const-> ActivityModelFnResultConst
// {
//     return {
//         rows(ln_activity_coefficients, ispecies, nspecies),
//         rows(ln_activities, ispecies, nspecies),
//         rows(partial_molar_volumes, ispecies, nspecies),
//         row(phase_molar_volumes, iphase, ispecies, nspecies),
//         row(phase_residual_molar_gibbs_energies, iphase, ispecies, nspecies),
//         row(phase_residual_molar_enthalpies, iphase, ispecies, nspecies),
//         row(phase_residual_molar_heat_capacities_cp, iphase, ispecies, nspecies),
//         row(phase_residual_molar_heat_capacities_cv, iphase, ispecies, nspecies)
//     };
// }

// } // namespace Reaktoro
