// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// The struct that stores the sensitivity of the optimal state with respect to a parameter
struct OptimumSensitivity
{
    /// The sensitivity of the primal variables `x` with respect to parameter `p`
    Vector dxdp;

    /// The sensitivity of the dual variables `y` with respect to parameter `p`
    Vector dydp;

    /// The sensitivity of the dual variables `z` with respect to parameter `p`
    Vector dzdp;
};

} // namespace Reaktoro
