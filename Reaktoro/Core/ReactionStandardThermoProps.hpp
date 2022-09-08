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
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// The primary standard thermodynamic properties of a chemical reaction.
/// In this type, the primary standard thermodynamic properties of a chemical
/// reaction are stored. This is a type to be used as the return type of a
/// function that calculates such properties of reactions. See
/// @ref ReactionStandardThermoModel. See also FormationReaction. Note there is no standard
/// molar volume property stored here. This is because a standard molar volume
/// model needs to be assigned to each individual Species object.
struct ReactionStandardThermoProps
{
    /// The standard molar Gibbs energy change @f$\Delta G^{\circ}@f$ of the reaction (in J/mol).
    real dG0;

    /// The standard molar enthalpy change @f$\Delta H^{\circ}@f$ of the reaction (in J/mol).
    real dH0;

    /// The standard molar isobaric heat capacity change @f$\Delta C_{P}^{\circ}@f$ of the reaction (in J/(mol·K)).
    real dCp0;
};

} // namespace Reaktoro
