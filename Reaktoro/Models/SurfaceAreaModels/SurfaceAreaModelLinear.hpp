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

/// Return a surface area model based on linear variation from a given molar surface area.
/// @param phase The name of the phase (as present in the chemical system) that forms the surface.
/// @param Abar The molar area of the surface (in m2/mol)
auto SurfaceAreaModelLinearMolar(String const& phase, real const& Abar) -> SurfaceAreaModel;

/// Return a surface area model based on linear variation from a given specific surface area.
/// @param phase The name of the phase (as present in the chemical system) that forms the surface.
/// @param Abar The specific area of the surface (in m2/kg)
auto SurfaceAreaModelLinearSpecific(String const& phase, real const& Abar) -> SurfaceAreaModel;

/// Return a surface area model based on linear variation from a given volumetric surface area.
/// @param phase The name of the phase (as present in the chemical system) that forms the surface.
/// @param Abar The volumetric area of the surface (in m2/m3)
auto SurfaceAreaModelLinearVolumetric(String const& phase, real const& Abar) -> SurfaceAreaModel;

/// Return a surface area model based on linear variation from a given normalized surface area.
/// This surface area model computes the surface area @eq{A} of a phase using
/// the equation:
///
/// @eqc{A=q\bar{A},}
///
/// where @eq{q} is defined as:
///
/// @eqc{q=\begin{cases}\text{phase amount} & \text{if \ensuremath{\hat{A}} normalized by amount}\\\text{phase mass} & \text{if \ensuremath{\hat{A}} normalized by mass}\\\text{phase volume} & \text{if \ensuremath{\hat{A}} normalized by volume}\end{cases}.}
///
/// @param phase The name of the phase (as present in the chemical system) that forms the surface.
/// @param Abar The normalized area of the surface.
/// @param unitAbar The unit of the normalized area of the surface (must be convertible to m2/mol, m2/kg, or m2/m3).
/// @ingroup Models
auto SurfaceAreaModelLinear(String const& phase, real Abar, Chars unitAbar) -> SurfaceAreaModel;

} // namespace Reaktoro
