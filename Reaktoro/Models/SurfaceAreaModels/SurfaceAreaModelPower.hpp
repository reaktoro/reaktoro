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
#include <Reaktoro/Core/SurfaceAreaModel.hpp>

namespace Reaktoro {

/// Return a surface area model for a phase based on power law that depends on a given initial surface area and phase amount.
/// This method is commonly used for solids and minerals (pure or solid
/// solutions). It computes the surface area @eq{A} of a phase with amount
/// @eq{q} using the equation:
///
/// @eqc{A(q)=A_{0}\left(\dfrac{q}{q_{0}}\right)^{p},}
///
/// where @eq{A_0} and @eq{q_0} are the initial surface area and amount of the
/// phase, and @eq{p} is a parameter. For uniformly dissolving spheres and
/// cubes, @eq{p = 2/3} (see @cite{Palandri2004}). For linearly varying surface
/// area, @eq{p = 1}.
/// @param phase The name of the phase (as present in the chemical system) that forms the surface.
/// @param A0 The initial surface area of the phase (in m2).
/// @param q0 The initial amount of the phase (in mol).
/// @param p The power parameter in the model.
/// @ingroup Models
auto SurfaceAreaModelPowerMolar(String const& phase, Param const& A0, Param const& q0, Param const& p) -> SurfaceAreaModel;

/// Return a surface area model for a phase based on power law that depends on a given initial surface area and phase mass.
/// This method is commonly used for solids and minerals (pure or solid
/// solutions). It computes the surface area @eq{A} of a phase with mass
/// @eq{q} using the equation:
///
/// @eqc{A(q)=A_{0}\left(\dfrac{q}{q_{0}}\right)^{p},}
///
/// where @eq{A_0} and @eq{q_0} are the initial surface area and mass of the
/// phase, and @eq{p} is a parameter. For uniformly dissolving spheres and
/// cubes, @eq{p = 2/3} (see @cite{Palandri2004}). For linearly varying surface
/// area, @eq{p = 1}.
/// @param phase The name of the phase (as present in the chemical system) that forms the surface.
/// @param A0 The initial surface area of the phase (in m2).
/// @param q0 The initial mass of the phase (in kg).
/// @param p The power parameter in the model.
/// @ingroup Models
auto SurfaceAreaModelPowerSpecific(String const& phase, Param const& A0, Param const& q0, Param const& p) -> SurfaceAreaModel;

/// Return a surface area model for a phase based on power law that depends on a given initial surface area and phase volume.
/// This method is commonly used for solids and minerals (pure or solid
/// solutions). It computes the surface area @eq{A} of a phase with volume
/// @eq{q} using the equation:
///
/// @eqc{A(q)=A_{0}\left(\dfrac{q}{q_{0}}\right)^{p},}
///
/// where @eq{A_0} and @eq{q_0} are the initial surface area and volume of the
/// phase, and @eq{p} is a parameter. For uniformly dissolving spheres and
/// cubes, @eq{p = 2/3} (see @cite{Palandri2004}). For linearly varying surface
/// area, @eq{p = 1}.
/// @param phase The name of the phase (as present in the chemical system) that forms the surface.
/// @param A0 The initial surface area of the phase (in m2).
/// @param q0 The initial volume of the phase (in m3).
/// @param p The power parameter in the model.
/// @ingroup Models
auto SurfaceAreaModelPowerVolumetric(String const& phase, Param const& A0, Param const& q0, Param const& p) -> SurfaceAreaModel;

/// Return a surface area model for a phase based on power law that depends on a given initial surface area and phase amount, mass or volume.
/// This method is identical to either one of below:
///   * SurfaceAreaModelPowerMolar (when `unitq0` is convertible to mol)
///   * SurfaceAreaModelPowerSpecific (when `unitq0` is convertible to kg)
///   * SurfaceAreaModelPowerVolumetric (when `unitq0` is convertible to m3)
/// @param phase The name of the phase (as present in the chemical system) that forms the surface.
/// @param A0 The initial surface area of the phase.
/// @param unitA0 The unit of the initial surface area (must be convertible to m2).
/// @param q0 The initial quantity of the phase.
/// @param unitq0 The unit of the initial quantity (must be convertible to `mol`, `kg`, or `m3`).
/// @param p The power parameter in the model.
/// @ingroup Models
auto SurfaceAreaModelPower(String const& mineral, real A0, Chars unitA0, real m0, Chars unitm0, Param const& q) -> SurfaceAreaModel;

} // namespace Reaktoro
