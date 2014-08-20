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
#include <ostream>
#include <tuple>

namespace Reaktor {

/**
 * Struct that defines a Cartesian point (x, y)
 */
struct Point
{
	/// The x-coordinate of the point
	double x;

	/// The y-coordinate of the point
	double y;

	/**
	 * Constructs a default @ref Point instance
	 */
	Point();

	/**
	 * Constructs a @ref Point instance
	 *
	 * @param x The x-coordinate of the point
	 * @param y The y-coordinate of the point
	 */
	Point(double x, double y);
};

/**
 * Overload operator that adds the points @c a and @c b
 *
 * @param a The left point in the addition
 * @param b The right point in the addition
 *
 * @return The addition of the two points @c a and @c b
 */
auto operator+(const Point& a, const Point& b) -> Point;

/**
 * Overload operator that subtracts the points @c a and @c b
 *
 * @param a The left point in the subtraction
 * @param b The right point in the subtraction
 *
 * @return The subtraction of the two points @c a and @c b
 */
auto operator-(const Point& a, const Point& b) -> Point;

/**
 * Overload operator that multiplies the point @c a by a scalar
 *
 * @param scalar The scalar used for the multiplication
 * @param a The point to be multiplied
 *
 * @return The multiplication of point @c a by @c scalar
 */
auto operator*(double scalar, const Point& a) -> Point;

/**
 * Overload operator that multiplies the point @c a by a scalar
 *
 * @param scalar The scalar used for the multiplication
 * @param a The point to be multiplied
 *
 * @return The multiplication of point @c a by @c scalar
 */
auto operator*(const Point& a, double scalar) -> Point;

/**
 * Overload operator that divides the point @c a by a scalar
 *
 * @param scalar The scalar used for the division
 * @param a The point to be divided
 *
 * @return The division of point @c a by @c scalar
 */
auto operator/(const Point& a, double scalar) -> Point;

/**
 * Calculates the squared distance between two Cartesian points
 *
 * @param a The Cartesian point @c a
 * @param b The Cartesian point @c b
 *
 * @return The squared distance between the Cartesian points @c a and @c b
 */
auto distanceSquared(const Point& a, const Point& b) -> double;

/**
 * Calculates the distance between two Cartesian points
 *
 * @param a The Cartesian point @c a
 * @param b The Cartesian point @c b
 *
 * @return The distance between the Cartesian points @c a and @c b
 */
auto distance(const Point& a, const Point& b) -> double;

/**
 * Outputs a Cartesian point (x, y)
 *
 * @param out The output stream
 * @param p The Cartesian point to be output
 *
 * @return The updated output stream
 */
auto operator<<(std::ostream& out, const Point& p) -> std::ostream&;

}  /* namespace Reaktor */
