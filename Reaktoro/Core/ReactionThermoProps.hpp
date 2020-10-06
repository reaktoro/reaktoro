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
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/Params.hpp>

namespace Reaktoro {

/// The primary standard thermodynamic properties of a chemical reaction.
/// In this type, the primary standard thermodynamic properties of a chemical
/// reaction are stored. This is a type to be used as the return type of a
/// function that calculates such properties of reactions. See @ref
/// ReactionThermoModel. See also @ref FormationReaction. Note there is no
/// standard molar volume property stored here. This is because a standard
/// molar volume model needs to be assigned to each individual @ref Species
/// object.
struct ReactionThermoProps
{
    /// The standard molar Gibbs energy change @f$\Delta G^{\circ}@f$ of the reaction (in J/mol).
    real dG0 = {};

    /// The standard molar enthalpy change @f$\Delta H^{\circ}@f$ of the reaction (in J/mol).
    real dH0 = {};
};

/// The arguments in a ReactionThermoModel function object.
struct ReactionThermoArgs
{
    /// The temperature for the calculation (in K)
    const real& T;

    /// The pressure for the calculation (in Pa)
    const real& P;

    /// The standard molar volume change @f$\Delta V^{\circ}@f$ of the reaction (in J/mol).
    const real& dV0;
};

/// Convenience macro to expand the structured bindings of a ReactionThermoArgs members.
#define ReactionThermoArgsDecl(args) \
    [[maybe_unused]] const auto& [T, P, dV0] = args;

/// The function type for calculation of standard thermodynamic properties of a reaction.
using ReactionThermoModel = Model<ReactionThermoProps(ReactionThermoArgs)>;

} // namespace Reaktoro
