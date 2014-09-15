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
#include <iostream>

namespace Reaktor {

// Forward declarations
struct WaterThermo;

struct WaterElectro
{
	WaterElectro();

	/// The dielectric constant of water
	double epsilon;

	/// The first-order partial derivative of the dielectric constant with respect to temperature
	double epsilonT;

	/// The first-order partial derivative of the dielectric constant with respect to pressure
	double epsilonP;

	/// The second-order partial derivative of the dielectric constant with respect to temperature
	double epsilonTT;

	/// The second-order partial derivative of the dielectric constant with respect to temperature and pressure
	double epsilonTP;

	/// The second-order partial derivative of the dielectric constant with respect to pressure
	double epsilonPP;

	/// The Born function \f$ Z\equiv-\frac{1}{\epsilon} \f$ (see Helgeson and Kirkham, 1974)
	double bornZ;

	/// The Born function \f$ Y\equiv\left[\frac{\partial Z}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	double bornY;

	/// The Born function \f$ Q\equiv\left[\frac{\partial Z}{\partial P}\right]_{T} \f$ (see Helgeson and Kirkham, 1974)
	double bornQ;

	/// The Born function \f$ N\equiv\left[\frac{\partial Q}{\partial P}\right]_{T} \f$ (see Helgeson and Kirkham, 1974)
	double bornN;

	/// The Born function \f$ U\equiv\left[\frac{\partial Q}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	double bornU;

	/// The Born function \f$ X\equiv\left[\frac{\partial Y}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	double bornX;
};

/// Output the electrostatic state of water
auto operator<<(std::ostream& out, const WaterElectro& we) -> std::ostream&;

// Calculate the electrostatic state of water using the model of Johnson and Norton (1991)
auto waterElectroState(double T, double P, const WaterThermo& wt) -> WaterElectro;

} // namespace Reaktor
