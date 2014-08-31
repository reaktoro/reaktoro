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

#include "Roots.hpp"

namespace Reaktor {

auto cubicRoots(double a, double b, double c, double d) -> CubicRoots
{
	const double xn = -b/(3*a);
	const double yn = a*xn*xn*xn + b*xn*xn + c*xn + d;

	const double d2    = (b*b - 3*a*c)/(9*a*a);
	const double h2    = 4*a*a*d2*d2*d2;
	const double Delta = yn*yn - h2;

	std::complex<double> x1, x2, x3;

	if(Delta > 0.0)
	{
		const double sqrtDelta = std::sqrt(Delta);

		const double operand1 = (-yn + sqrtDelta)/(2*a);
		const double operand2 = (-yn - sqrtDelta)/(2*a);

		const double alpha = xn +
			pow(std::abs(operand1), 1.0/3) * std::abs(operand1)/operand1 +
			pow(std::abs(operand2), 1.0/3) * std::abs(operand2)/operand2;

		const double discr = b*b - 4*a*c - 2*a*b*alpha - 3*pow(a*alpha, 2);
		const double aux   = std::sqrt(-discr);

		x1 = {alpha, 0.0};
		x2 = {(-b - a*alpha)/(2*a), -aux/(2*a)};
		x3 = {(-b - a*alpha)/(2*a),  aux/(2*a)};
	}
	else
	{
		const double pi    = 3.14159265359;
		const double delta = std::sqrt((b*b - 3*a*c)/(9*a*a));
		const double h     = 2*a*delta*delta*delta;
		const double theta = acos(-yn/h)/3;

		const double alpha = xn + 2*delta*cos(theta);
		const double beta  = xn + 2*delta*cos(2*pi/3 - theta);
		const double gamma = xn + 2*delta*cos(2*pi/3 + theta);

		x1 = {alpha, 0.0};
		x2 = {beta,  0.0};
		x3 = {gamma, 0.0};
	}

	return std::make_tuple(x1, x2, x3);
}

} // namespace Reaktor
