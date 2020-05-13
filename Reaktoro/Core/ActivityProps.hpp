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
#include <any>

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// The base type for the primary activity and excess thermodynamic properties of a phase.
/// @see ActivityPropsFn, ActivityArgs
template<typename Real, typename Array>
struct ActivityPropsBase
{
    /// The excess molar volume of the phase (in m3/mol).
    Real Vex;

    /// The temperature derivative of the excess molar volume at constant pressure (in m3/(mol*K)).
    Real VexT;

    /// The pressure derivative of the excess molar volume at constant temperature (in m3/(mol*Pa)).
    Real VexP;

    /// The excess molar Gibbs energy of the phase (in units of J/mol).
    Real Gex;

    /// The excess molar enthalpy of the phase (in units of J/mol).
    Real Hex;

    /// The excess molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    Real Cpex;

    /// The excess molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    Real Cvex;

    /// The activity coefficients (natural log) of the species in the phase.
    Array ln_g;

    /// The activities (natural log) of the species in the phase.
    Array ln_a;

    /// Assign a common value to all properties in this ActivityPropsBase object.
    auto operator=(real value) -> ActivityPropsBase&
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

    /// Convert this ActivityPropsBase object into another.
    template<typename RX, typename AX>
    operator ActivityPropsBase<RX, AX>()
    {
        return { Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a };
    }

    /// Convert this ActivityPropsBase object into another.
    template<typename RX, typename AX>
    operator ActivityPropsBase<RX, AX>() const
    {
        return { Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a };
    }

    /// Create a ActivityPropsBase object with given number of species.
    /// This static method is needed, instead of a constructor, which would
    /// prevent aggregate initialization of this struct.
    static auto create(Index numspecies) -> ActivityPropsBase<Real, Array>
    {
        ActivityPropsBase<Real, Array> props = {};
        props.Vex  = 0.0;
        props.VexT = 0.0;
        props.VexP = 0.0;
        props.Gex  = 0.0;
        props.Hex  = 0.0;
        props.Cpex = 0.0;
        props.Cvex = 0.0;
        props.ln_g = Array::Zero(numspecies);
        props.ln_a = Array::Zero(numspecies);
        return props;
    }
};

/// The activity and excess thermodynamic properties of a phase.
using ActivityProps = ActivityPropsBase<real, ArrayXr>;

/// The non-const view to the activity and excess thermodynamic properties of a phase.
using ActivityPropsRef = ActivityPropsBase<real&, ArrayXrRef>;

/// The const view to the activity and excess thermodynamic properties of a phase.
using ActivityPropsConstRef = ActivityPropsBase<const real&, ArrayXrConstRef>;

/// The arguments in an function for calculation of activity properties of a phase.
/// @see ActivityPropsFn, ActivityProps
struct ActivityArgs
{
    /// The temperature for the calculation (in K).
    const real& T;

    /// The pressure for the calculation (in Pa).
    const real& P;

    /// The mole fractions of the species in the phase.
    ArrayXrConstRef x;

    /// The extra arguments for the activity model evaluation whose type is only known at runtime.
    Vec<Any>& extra;
};

/// The function type for the calculation of activity and excess thermodynamic properties of a phase.
using ActivityPropsFn = Fn<void(ActivityPropsRef, ActivityArgs)>;

} // namespace Reaktoro
