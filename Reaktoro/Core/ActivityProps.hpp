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

// C++ includes
#include <any>

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/Model.hpp>

namespace Reaktoro {

/// The base type for the primary activity and excess thermodynamic properties of a phase.
/// @see ActivityModel, ActivityArgs
template<typename Real, typename Array, typename Extra>
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

    /// The extra data produced by an activity model that may be reused by subsequent models within a chained activity model.
    Extra extra;

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
    template<typename RX, typename AX, typename EX>
    auto operator=(const ActivityPropsBase<RX, AX, EX>& other) -> ActivityPropsBase&
    {
        Vex   = other.Vex;
        VexT  = other.VexT;
        VexP  = other.VexP;
        Gex   = other.Gex;
        Hex   = other.Hex;
        Cpex  = other.Cpex;
        Cvex  = other.Cvex;
        ln_g  = other.ln_g;
        ln_a  = other.ln_a;
        extra = other.extra;
        return *this;
    }

    /// Convert this ActivityPropsBase object into another.
    template<typename RX, typename AX, typename EX>
    operator ActivityPropsBase<RX, AX, EX>()
    {
        return { Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a, extra };
    }

    /// Convert this ActivityPropsBase object into another.
    template<typename RX, typename AX, typename EX>
    operator ActivityPropsBase<RX, AX, EX>() const
    {
        return { Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a, extra };
    }

    /// Create a ActivityPropsBase object with given number of species.
    /// This static method is needed, instead of a constructor, which would
    /// prevent aggregate initialization of this struct.
    static auto create(Index numspecies) -> ActivityPropsBase<Real, Array, Extra>
    {
        ActivityPropsBase<Real, Array, Extra> props = {};
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
using ActivityProps = ActivityPropsBase<real, ArrayXr, Vec<Any>>;

/// The non-const view to the activity and excess thermodynamic properties of a phase.
using ActivityPropsRef = ActivityPropsBase<real&, ArrayXrRef, Vec<Any>&>;

/// The const view to the activity and excess thermodynamic properties of a phase.
using ActivityPropsConstRef = ActivityPropsBase<const real&, ArrayXrConstRef, const Vec<Any>&>;

} // namespace Reaktoro
