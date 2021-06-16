// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Core/StandardThermoProps.hpp>

namespace Reaktoro {

/// The parameters in the constant model for calculating standard thermodynamic properties of species.
struct StandardThermoModelParamsConstant
{
    /// The constant standard molar Gibbs energy @f$G^{\circ}@f$ of the species (in J/mol).
    Param G0 = Param("G0", 0.0);

    /// The constant standard molar enthalpy @f$H^{\circ}@f$ of the species (in J/mol).
    Param H0 = Param("H0", 0.0);

    /// The constant standard molar volume @f$V^{\circ}@f$ of the species (in m3/mol).
    Param V0 = Param("V0", 0.0);

    /// The constant standard molar isobaric heat capacity @f$C_{P}^{\circ}@f$ of the species (in J/(mol·K)).
    Param Cp0 = Param("Cp0", 0.0);

    /// The constant standard molar isochoric heat capacity @f$C_{V}^{\circ}@f$ of the species (in J/(mol·K)).
    Param Cv0 = Param("Cv0", 0.0);
};

/// Return a function that calculates thermodynamic properties of a species using a constant model for its standard properties.
auto StandardThermoModelConstant(const StandardThermoModelParamsConstant& params) -> StandardThermoModel;

} // namespace Reaktoro
