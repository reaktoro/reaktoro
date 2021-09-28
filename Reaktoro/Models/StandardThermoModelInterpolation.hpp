// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/StandardThermoProps.hpp>

namespace Reaktoro {

/// The parameters in the Maier-Kelley model for calculating standard thermodynamic properties of fluid and solid species.
struct StandardThermoModelParamsInterpolation
{
    /// The temperatures at which known thermodynamic data is known (in K).
    Vec<double> temperatures;

    /// The pressures at which known thermodynamic data is known (in Pa).
    Vec<double> pressures;

    /// The standard molar Gibbs energies @f$G^{\circ}@f$ of formation of the species (in J/mol).
    Vec<Vec<double>> G0;

    /// The standard molar enthalpies @f$H^{\circ}@f$ of formation of the species (in J/mol).
    Vec<Vec<double>> H0;

    /// The standard molar volumes @f$V^{\circ}@f$ of the species (in m³/mol).
    Vec<Vec<double>> V0;

    /// The temperature derivative of the standard molar volumes @f$\partial V^{\circ}/\partial T@f$ of the species (in m³/(mol·K)).
    Vec<Vec<double>> VT0;

    /// The pressure derivative of the standard molar volumes @f$\partial V^{\circ}/\partial P@f$ of the species (in m³/(mol·K)).
    Vec<Vec<double>> VP0;

    /// The standard molar isobaric heat capacities @f$C_{P}^{\circ}@f$ of the species (in J/(mol·K)).
    Vec<Vec<double>> Cp0;

    /// The reference pressure used for volume correction of @f$G^{\circ}@f$ (in Pa).
    double Pref = 1e5;
};

/// Return a function that calculates thermodynamic properties of a species using the Maier-Kelley model.
auto StandardThermoModelInterpolation(const StandardThermoModelParamsInterpolation& params) -> StandardThermoModel;

} // namespace Reaktoro
