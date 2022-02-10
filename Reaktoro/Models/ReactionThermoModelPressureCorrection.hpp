// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Core/ReactionThermoProps.hpp>

namespace Reaktoro {

/// Return a function that calculates pressure correction for standard Gibbs energy and enthalpy a reaction.
///
/// In this model, the standard Gibbs energy, enthalpy, and isobatic heat capacity of reaction are computed using:
///
/// @eqc{\Delta G^{\circ}=\Delta G_{\mathrm{base}}^{\circ}+\Delta V^{\circ}(P-P_{r})}
/// @eqc{\Delta H^{\circ}=\Delta H_{\mathrm{base}}^{\circ}+\Delta V^{\circ}(P-P_{r})}
/// @eqc{\Delta C_P^{\circ}=\Delta C_{P,\mathrm{base}}^{\circ}}
///
/// where @eq{\Delta G_{\mathrm{base}}^{\circ}}, @eq{\Delta H_{\mathrm{base}}^{\circ}},
/// and @eq{\Delta C_{P,\mathrm{base}}^{\circ}} denote standard properties computed using
/// a given base thermodynamic model function of the reaction.
/// @param Pr The reference pressure for the pressure correction (in Pa).
auto ReactionThermoModelPressureCorrection(Param Pr) -> ReactionThermoModel;

} // namespace Reaktoro
