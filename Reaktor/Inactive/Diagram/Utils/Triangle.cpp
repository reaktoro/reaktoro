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

#include "Triangle.hpp"

// C++ includes
#include <cmath>

namespace Reaktor {

Triangle::Triangle()
{
}

Triangle::Triangle(const Point& a, const Point& b, const Point& c)
: a(a), b(b), c(c)
{
}

auto length(const Triangle& t) -> double
{
	const double dist_ab = distance(t.a, t.b);
	const double dist_bc = distance(t.b, t.c);
	const double dist_ca = distance(t.c, t.a);

	return std::min(dist_ab, std::min(dist_bc, dist_ca));
}

}  /* namespace Reaktor */
