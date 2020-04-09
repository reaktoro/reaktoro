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

// C++ includes
#include <cmath>
#include <functional>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// The activity and excess thermodynamic properties of a phase.
/// @see ActivityPropsFn, StandardThermoPropsFn, StandardThermoProps
struct ActivityProps
{
    /// The excess molar volume of the phase (in m3/mol).
    real& Vex;

    /// The temperature derivative of the excess molar volume at constant pressure (in m3/(mol*K)).
    real& VexT;

    /// The pressure derivative of the excess molar volume at constant temperature (in m3/(mol*Pa)).
    real& VexP;

    /// The excess molar Gibbs energy of the phase (in units of J/mol).
    real& Gex;

    /// The excess molar enthalpy of the phase (in units of J/mol).
    real& Hex;

    /// The excess molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    real& Cpex;

    /// The excess molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    real& Cvex;

    /// The activity coefficients (natural log) of the species in the phase.
    ArrayXrRef ln_g;

    /// The activities (natural log) of the species in the phase.
    ArrayXrRef ln_a;
};

/// The function type for calculation of activity and excess thermodynamic properties of a phase.
using ActivityPropsFn = std::function<void(ActivityProps, real, real, ArrayXrConstRef)>;

/// The function type for activity model of a phase.
using ActivityModel = std::function<ActivityPropsFn(const SpeciesList&)>;

} // namespace Reaktoro
