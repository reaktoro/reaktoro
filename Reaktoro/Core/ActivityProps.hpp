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

// C++ includes
#include <any>

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/TypeOp.hpp>
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>

namespace Reaktoro {

/// The base type for the primary activity and corrective thermodynamic
/// properties of a phase. Thermodynamic properties for a phase, such as
/// internal energy, enthalpy, Gibbs energy, entropy, and volume can be broken
/// down into *ideal* and *corrective* contributions. Let us denote by @eq{M}
/// one of these properties. The previous statement implies that:
///
/// @eqc{M=\sum_{i=1}^{\mathrm{N}}x_{i}M_{i}^{\circ}+M^{\mathrm{x}},}
///
/// with the first term being the *ideal contribution* and the second term,
/// @eq{M^{\mathrm{x}}}, being the *corrective contribution*, where @eq{x_{i}}
/// and @eq{M_{i}^{\circ}} are, respectively, the mole fraction and respective
/// standard molar property of the \eq{i}-th species in the phase.
///
/// Thus, corrective thermodynamic properties are those that need to be added
/// to their ideal counterpart to obtain more correct values that take into
/// account the non-ideal thermodynamic behavior of phases. The word
/// *corrective* is adopted in %Reaktoro to mean either *excess* or *residual*
/// properties.
///
/// @note The corrective property \eq{M^{\mathrm{x}}} may sometimes be the
/// *actual complete property* of the phase, i.e., \eq{M^{\mathrm{x}} \equiv
/// M}. For example, for gaseous phases in which the partial molar volumes of
/// the species are conventionally zero, @eq{V_{i}^{\circ}=0}, the corrective
/// molar volume, @eq{V^{\mathrm{x}}}, must be set as the total molar volume of
/// the phase. @see ActivityModel, ActivityArgs
template<template<typename> typename TypeOp>
struct ActivityPropsBase
{
    /// The corrective molar volume of the phase (in m3/mol).
    TypeOp<real> Vx;

    /// The temperature derivative of the corrective molar volume at constant pressure (in m3/(mol*K)).
    TypeOp<real> VxT;

    /// The pressure derivative of the corrective molar volume at constant temperature (in m3/(mol*Pa)).
    TypeOp<real> VxP;

    /// The corrective molar Gibbs energy of the phase (in units of J/mol).
    TypeOp<real> Gx;

    /// The corrective molar enthalpy of the phase (in units of J/mol).
    TypeOp<real> Hx;

    /// The corrective molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    TypeOp<real> Cpx;

    /// The activity coefficients (natural log) of the species in the phase.
    TypeOp<ArrayXr> ln_g;

    /// The activities (natural log) of the species in the phase.
    TypeOp<ArrayXr> ln_a;

    /// The state of matter of the phase.
    TypeOp<StateOfMatter> som;

    /// The extra data produced by an activity model that may be reused by subsequent models within a chained activity model.
    TypeOp<Map<String, Any>> extra;

    /// Assign a common value to all properties in this ActivityPropsBase object.
    auto operator=(real value) -> ActivityPropsBase&
    {
        Vx   = value;
        VxT  = value;
        VxP  = value;
        Gx   = value;
        Hx   = value;
        Cpx  = value;
        ln_g = value;
        ln_a = value;
        return *this;
    }

    /// Convert this ActivityPropsBase object into another.
    template<template<typename> typename OtherTypeOp>
    auto operator=(const ActivityPropsBase<OtherTypeOp>& other) -> ActivityPropsBase&
    {
        Vx    = other.Vx;
        VxT   = other.VxT;
        VxP   = other.VxP;
        Gx    = other.Gx;
        Hx    = other.Hx;
        Cpx   = other.Cpx;
        ln_g  = other.ln_g;
        ln_a  = other.ln_a;
        som   = other.som;
        extra = other.extra;
        return *this;
    }

    /// Convert this ActivityPropsBase object into another.
    template<template<typename> typename OtherTypeOp>
    operator ActivityPropsBase<OtherTypeOp>()
    {
        return { Vx, VxT, VxP, Gx, Hx, Cpx, ln_g, ln_a, som, extra };
    }

    /// Convert this ActivityPropsBase object into another.
    template<template<typename> typename OtherTypeOp>
    operator ActivityPropsBase<OtherTypeOp>() const
    {
        return { Vx, VxT, VxP, Gx, Hx, Cpx, ln_g, ln_a, som, extra };
    }

    /// Create a ActivityPropsBase object with given number of species.
    /// This static method is needed, instead of a constructor, which would
    /// prevent aggregate initialization of this struct.
    static auto create(Index numspecies) -> ActivityPropsBase<TypeOp>
    {
        ActivityPropsBase<TypeOp> props = {};
        props.Vx  = 0.0;
        props.VxT = 0.0;
        props.VxP = 0.0;
        props.Gx  = 0.0;
        props.Hx  = 0.0;
        props.Cpx = 0.0;
        props.ln_g = ArrayXr::Zero(numspecies);
        props.ln_a = ArrayXr::Zero(numspecies);
        return props;
    }
};

/// The activity and corrective thermodynamic properties of a phase.
using ActivityProps = ActivityPropsBase<TypeOpIdentity>;

/// The non-const view to the activity and corrective thermodynamic properties of a phase.
using ActivityPropsRef = ActivityPropsBase<TypeOpRef>;

/// The const view to the activity and corrective thermodynamic properties of a phase.
using ActivityPropsConstRef = ActivityPropsBase<TypeOpConstRef>;

} // namespace Reaktoro
