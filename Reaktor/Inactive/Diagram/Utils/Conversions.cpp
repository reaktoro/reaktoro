/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "Conversions.hpp"

// C++ includes
#include <cmath>

// Reaktor includes
#include <Reaktor/Diagram/Utils/BarycentricPoint.hpp>
#include <Reaktor/Diagram/Utils/Point.hpp>

namespace Reaktor {

auto barycentricPoint(const Point& p) -> BarycentricPoint
{
	const double x = p.x;
	const double y = p.y;

	const double sqrt3 = std::sqrt(3.0);

	const double B = x - y/sqrt3;
	const double C = 2 * y/sqrt3;
	const double A = 1.0 - B - C;

	return (0 <= B and B <= 1 and 0 <= C and C <= 1) ?
		BarycentricPoint(A, B, C) : BarycentricPoint();
}

auto point(const BarycentricPoint& q) -> Point
{
	const double B = q.B;
	const double C = q.C;

	const double sqrt3 = std::sqrt(3.0);

	const double x = B + 0.5 * C;
	const double y = sqrt3/2.0 * C;

	return Point(x, y);
}

}  /* namespace Reaktor */
