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

#include "BarycentricPoint.hpp"

// C++ includes
#include <cmath>

namespace Reaktor {

BarycentricPoint::BarycentricPoint()
: A(0), B(0), C(0)
{
}

BarycentricPoint::BarycentricPoint(double A, double B)
: A(A), B(B), C(1 - A - B)
{
}

BarycentricPoint::BarycentricPoint(double A, double B, double C)
: A(A), B(B), C(C)
{
}

auto operator<<(std::ostream& out, const BarycentricPoint& q) -> std::ostream&
{
    out << "(" << q.A << ", " << q.B << ", " << q.C << ")";

    return out;
}

}  /* namespace Reaktor */
