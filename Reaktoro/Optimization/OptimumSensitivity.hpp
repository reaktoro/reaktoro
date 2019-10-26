// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// A type where the sensitivity derivatives of the computed optimum state is stored.
struct OptimumSensitivity
{
    /// The sensitivity derivatives of the primal solution *x* with respect to parameters *p*.
    Matrix dxdp;

    /// The sensitivity derivatives of the dual solution *y* with respect to parameters *p*.
    Matrix dydp;

    /// The sensitivity derivatives of the dual solution *z* with respect to parameters *p*.
    Matrix dzdp;

    /// The sensitivity derivatives of the dual solution *w* with respect to parameters *p*.
    Matrix dwdp;
};

} // namespace Reaktoro
