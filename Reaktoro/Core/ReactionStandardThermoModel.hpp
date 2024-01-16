// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/ReactionStandardThermoProps.hpp>

namespace Reaktoro {

/// The arguments in a ReactionStandardThermoModel function object.
struct ReactionStandardThermoModelArgs
{
    /// The temperature for the calculation (in K)
    real const& T;

    /// The pressure for the calculation (in Pa)
    real const& P;

    /// The standard molar volume change @f$\Delta V^{\circ}@f$ of the reaction (in J/mol).
    real const& dV0;
};

/// The function type for calculation of standard thermodynamic properties of a reaction.
using ReactionStandardThermoModel = Model<ReactionStandardThermoProps(ReactionStandardThermoModelArgs)>;

} // namespace Reaktoro
