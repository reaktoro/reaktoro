// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// C++ includes
#include <vector>

namespace Reaktoro {

class LagrangeInterpolator
{
public:
    LagrangeInterpolator();

    LagrangeInterpolator(const std::vector<double>& xi, const std::vector<double>& yi, unsigned order = 1);

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
