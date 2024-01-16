// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "SurfaceAreaModelConstant.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>

namespace Reaktoro {

auto SurfaceAreaModelConstant(real const& A0) -> SurfaceAreaModel
{
    return [=](ChemicalProps const& props) -> real { return A0; };
}

auto SurfaceAreaModelConstant(real const& A0, Chars unitA0) -> SurfaceAreaModel
{
    auto A0m2 = units::convert(A0, unitA0, "m2");
    return SurfaceAreaModelConstant(A0m2);
}

} // namespace Reaktoro
