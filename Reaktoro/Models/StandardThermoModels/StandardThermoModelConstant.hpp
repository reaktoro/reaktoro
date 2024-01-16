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
#include <Reaktoro/Core/StandardThermoModel.hpp>

namespace Reaktoro {

/// The parameters in the constant model for calculating standard thermodynamic properties of species.
struct StandardThermoModelParamsConstant
{
    /// The constant standard molar Gibbs energy @eq{G^{\circ}} of the species (in J/mol).
    real G0;

    /// The constant standard molar enthalpy @eq{H^{\circ}} of the species (in J/mol).
    real H0;

    /// The constant standard molar volume @eq{V^{\circ}} of the species (in m³/mol).
    real V0;

    /// The temperature derivative of the standard molar volume @eq{\partial V^{\circ} / \partial T} of the species (in m³/(mol·K)).
    real VT0;

    /// The pressure derivative of the standard molar volume @eq{\partial V^{\circ} / \partial P} of the species (in m³/(mol·K)).
    real VP0;

    /// The constant standard molar isobaric heat capacity @eq{C_{P}^{\circ}} of the species (in J/(mol·K)).
    real Cp0;
};

/// Return a function that calculates thermodynamic properties of a species using a constant model for its standard properties.
auto StandardThermoModelConstant(const StandardThermoModelParamsConstant& params) -> StandardThermoModel;

} // namespace Reaktoro
