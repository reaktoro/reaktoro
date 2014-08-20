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

#pragma once

// Reaktor includes
#include <Reaktor/Diagram/Utils/Point.hpp>

namespace Reaktor {

/**
 * Struct that defines a triangle in the Cartesian space
 */
struct Triangle
{
	/// The first point of the triangle
	Point a;

	/// The second point of the triangle
	Point b;

	/// The third point of the triangle
	Point c;

	/**
	 * Constructs a default @ref Triangle instance
	 */
	Triangle();

	/**
	 * Constructs a @ref Triangle instance
	 *
	 * @param a The first point of the triangle
	 * @param b The second point of the triangle
	 * @param c The third point of the triangle
	 */
	Triangle(const Point& a, const Point& b, const Point& c);
};

/**
 * Calculates the length of a triangle
 *
 * The length of a triangle is calculated as the minimum
 * length of its sides.
 *
 * @param t The triangle as an instance of @ref Triangle
 *
 * @return The length of the triangle @c t
 */
auto length(const Triangle& triangle) -> double;

}  /* namespace Reaktor */
