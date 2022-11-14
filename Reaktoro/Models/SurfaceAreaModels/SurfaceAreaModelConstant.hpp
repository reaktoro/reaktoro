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

/// Return a constant surface area model.
/// @param A0 The value of the constant surface area (in m2).
/// @ingroup Models
auto SurfaceAreaModelConstant(Param const& A0) -> SurfaceAreaModel;

/// Return a constant surface area model.
/// @param A0 The value of the constant surface area.
/// @param unitA0 The unit of the surface area value (must be convertible to m2).
/// @ingroup Models
auto SurfaceAreaModelConstant(real const& A0, Chars unitA0) -> SurfaceAreaModel;

} // namespace Reaktoro
