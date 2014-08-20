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

namespace Reaktor {

/**
 * Struct that defines a barycentric point (A, B, C)
 *
 * The barycentric point (A, B, C) has the property of
 * 0 <= A, B, C <= 1 and A + B + C = 1.
 */
struct BarycentricPoint
{
	/// The A-coordinate of the barycentric point
	double A;

	/// The B-coordinate of the barycentric point
	double B;

	/// The C-coordinate of the barycentric point
	double C;

	/**
	 * Constructs a default @ref BarycentricPoint instance
	 */
	BarycentricPoint();

	/**
	 *  Constructs a @ref BarycentricPoint instance
	 *
	 *	The C-coordinate of the barycentric point is
	 *	automatically calculated as C = 1 - A - B.
	 *
	 *	@param A The A-coordinate of the barycentric point
	 *	@param B The B-coordinate of the barycentric point
	 */
	BarycentricPoint(double A, double B);

	/**
	 *  Constructs a @ref BarycentricPoint instance
	 *
	 *	@param A The A-coordinate of the barycentric point
	 *	@param B The B-coordinate of the barycentric point
	 *	@param C The C-coordinate of the barycentric point
	 */
	BarycentricPoint(double A, double B, double C);
};

/**
 * Outputs a barycentric point (A, B, C)
 *
 * @param out The output stream
 * @param q The barycentric point to be output
 *
 * @return The updated output stream
 */
auto operator<<(std::ostream& out, const BarycentricPoint& q) -> std::ostream&;

}  /* namespace Reaktor */
