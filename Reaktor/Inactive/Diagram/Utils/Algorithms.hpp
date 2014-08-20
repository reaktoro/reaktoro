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

// C++ includes
#include <vector>

namespace Reaktor {

// Reaktor forward declarations
struct Segment;
struct Point;

/**
 * Checks if one segment cross another
 *
 * @param s1 The first segment as an instance of @ref Segment
 * @param s2 The second segment as an instance of @ref Segment
 */
auto cross(const Segment& s1, const Segment& s2) -> bool;

/**
 * Checks if a segment cross a path
 *
 * @param s The segment as an instance of @ref Segment
 * @param path The set of ordered points describing a path
 */
auto cross(const Segment& s, const std::vector<Point>& path) -> bool;

/**
 * Sorts a set of points in order of distance from a reference point
 *
 * @param points The points to be sorted, as @ref Point instances
 * @param ref The reference point from which the distance is calculated
 *
 * @return The sorted points in order of distance from @c ref
 */
auto sort(const std::vector<Point>& points, const Point& ref) -> std::vector<Point>;

/**
 * Finds a continuous path from an set of unordered points
 *
 * @param points The set of points from which the path is found
 *
 * @return A ordered set of points describing a path
 */
auto findPath(const std::vector<Point>& points) -> std::vector<Point>;

}  /* namespace Reaktor */
