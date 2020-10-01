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
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// The primary standard thermodynamic properties of a chemical reaction.
/// In this type, the primary standard thermodynamic properties of a chemical
/// reaction are stored. This is a type to be used as the return type of a
/// function that calculates such properties of reactions. See @ref
/// ReactionThermoPropsFn. See also @ref FormationReaction. Note there is no
/// standard molar volume property stored here. This is because a standard
/// molar volume model needs to be assigned to each individual @ref Species
/// object.
struct ReactionThermoProps
{
    /// The standard molar Gibbs energy @f$\Delta G^{\circ}@f$ of the reaction (in J/mol).
    real dG0 = {};

    /// The standard molar enthalpy @f$\Delta H^{\circ}@f$ of the reaction (in J/mol).
    real dH0 = {};
};

/// The function type for calculation of standard thermodynamic properties of a reaction.
/// @param T The temperature for the calculation (in K)
/// @param P The pressure for the calculation (in Pa)
/// @return The standard thermodynamic properties of the reaction
using ReactionThermoPropsFn = Fn<ReactionThermoProps(real T, real P)>;

} // namespace Reaktoro
