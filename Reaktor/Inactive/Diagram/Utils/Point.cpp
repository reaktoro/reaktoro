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

#include "Point.hpp"

// C++ includes
#include <cmath>

namespace Reaktor {

Point::Point()
: x(0), y(0)
{
}

Point::Point(double x, double y)
: x(x), y(y)
{
}

auto operator+(const Point& a, const Point& b) -> Point
{
	return Point(a.x + b.x, a.y + b.y);
}

auto operator-(const Point& a, const Point& b) -> Point
{
	return Point(a.x - b.x, a.y - b.y);
}

auto operator*(double scalar, const Point& a) -> Point
{
	return Point(scalar * a.x, scalar * a.y);
}

auto operator*(const Point& a, double scalar) -> Point
{
	return scalar * a;
}

auto operator/(const Point& a, double scalar) -> Point
{
	return 1.0/scalar * a;
}

auto distanceSquared(const Point& a, const Point& b) -> double
{
	return std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2);
}

auto distance(const Point& a, const Point& b) -> double
{
	return std::sqrt(distanceSquared(a, b));
}

auto operator<<(std::ostream& out, const Point& p) -> std::ostream&
{
    out << "(" << p.x << ", " << p.y << ")";

    return out;
}

}  /* namespace Reaktor */
