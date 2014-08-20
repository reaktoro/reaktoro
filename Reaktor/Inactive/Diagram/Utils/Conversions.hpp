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

namespace Reaktor {

// Reaktor forward declarations
struct BarycentricPoint;
struct Point;

/**
 * Converts a Cartesian point into a barycentric point
 *
 * If the given Cartesian point (x, y) falls outside the barycentric triangle,
 * the bull barycentric point (0, 0, 0) is returned.
 *
 * @param p A Cartesian point (x, y) as an instance of @ref Point
 *
 * @return A barycentric point (A, B, C) as an instance of @ref BarycentricPoint
 */
auto barycentricPoint(const Point& p) -> BarycentricPoint;

/**
 * Converts a barycentric point into a Cartesian point
 *
 * @param q A barycentric point (A, B, C) as an instance of @ref BarycentricPoint
 *
 * @return A Cartesian point (x, y) as an instance of @ref Point
 */
auto point(const BarycentricPoint& q) -> Point;

}  /* namespace Reaktor */
