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

#pragma once

// C++ includes
#include <vector>

namespace Reaktoro {

/// A class used to calculate interpolation of data in one dimension in any order.
class LagrangeInterpolator
{
public:
    /// Construct a default LagrangeInterpolator instance.
    LagrangeInterpolator();

    /// Construct a LagrangeInterpolator instance.
    /// @param xp The \eq{x}-coordinate points of the \eq{y}-data to be interpolated.
    /// @param yp The \eq{y}-data points to be interpolated.
    /// @param order The order of the interpolation.
    LagrangeInterpolator(const std::vector<double>& xp, const std::vector<double>& yp, unsigned order = 1);

    /// Return the interpolation of \eq{y} at a given coordinate \eq{x}.
    auto operator()(double x) const -> double;

private:
    /// The interpolation points on the x-axis
    std::vector<double> xi;

    /// The interpolation points on the y-axis
    std::vector<double> yi;

    /// The order of the interpolation
    unsigned order;
};

} // namespace Reaktoro
