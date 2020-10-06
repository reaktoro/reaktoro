// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

/// Return a function that calculates thermodynamic properties of a reaction using a constant lg(K) model.
///
/// In this model, the equilibrium constant of the reaction remains constant at all temperatures.
/// is computed at a given temperature with:
///
/// @eqc{\log_10 K=\log_10 K_{0},}
///
/// where @eq{\ln K_{0}} is a given equilibrium constant.
/// The standard Gibbs energy of reaction is computed using:
///
/// @eqc{\Delta G^{\circ}=-RT\ln K_{0}}
///
/// while the standard enthalpy of reaction is assumed zero:
///
/// @eqc{\Delta H^{\circ}=0.}
///
/// @param lgK0 The equilibrium constant (log base 10) of the reaction.
auto ReactionThermoModelConstLgK(real lgK0) -> ReactionThermoModel;

} // namespace Reaktoro
