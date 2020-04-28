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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// The activity and excess thermodynamic properties of a phase.
/// @see ActivityPropsFn, ActivityArgs
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

    /// Assign a common value to all properties in this ActivityProps object.
    auto operator=(real value) -> ActivityProps&
    {
        Vex  = value;
        VexT = value;
        VexP = value;
        Gex  = value;
        Hex  = value;
        Cpex = value;
        Cvex = value;
        ln_g = value;
        ln_a = value;
        return *this;
    }
};

/// The arguments in an function for calculation of activity properties of a phase.
/// @see ActivityPropsFn, ActivityProps
struct ActivityArgs
{
    /// The temperature for the calculation (in K).
    real const& T;

    /// The pressure for the calculation (in Pa).
    real const& P;

    /// The mole fractions of the species in the phase.
    ArrayXrConstRef x;

    /// The extra arguments for the activity model evaluation whose type is only known at runtime.
    Vec<std::any>& extra;
};

/// The function type for calculation of activity and excess thermodynamic properties of a phase.
using ActivityPropsFn = Fn<void(ActivityProps, ActivityArgs)>;

/// The function type for the creation of activity model function of a phase.
using ActivityModel = Fn<ActivityPropsFn(const SpeciesList&)>;

} // namespace Reaktoro
