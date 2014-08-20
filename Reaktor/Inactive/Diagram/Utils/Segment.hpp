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
 * Struct that defines a segment of line in the Cartesian space
 */
struct Segment
{
	/// The start point of the segment
	Point a;

	/// The end point of the segment
	Point b;

	/**
	 * Constructs a default @ref Segment instance
	 */
	Segment();

	/**
	 * Constructs a @ref Segment instance
	 *
	 * @param a The start point of the segment
	 * @param b The end point of the segment
	 */
	Segment(const Point& a, const Point& b);
};

/**
 * Calculates the length of a segment
 *
 * The length of a segment is calculated as the distance
 * between its two vertices.
 *
 * @param s The segment as an instance of @ref Segment
 *
 * @return The length of the segment @c s
 */
auto length(const Segment& segment) -> double;

}  /* namespace Reaktor */
