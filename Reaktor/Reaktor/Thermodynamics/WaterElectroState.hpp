// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

namespace Reaktor {

struct WaterElectroState
{
	/// The dielectric constant of water
	double epsilon = 0.0;

	/// The first-order partial derivative of the dielectric constant with respect to temperature
	double epsilonT = 0.0;

	/// The first-order partial derivative of the dielectric constant with respect to pressure
	double epsilonP = 0.0;

	/// The second-order partial derivative of the dielectric constant with respect to temperature
	double epsilonTT = 0.0;

	/// The second-order partial derivative of the dielectric constant with respect to temperature and pressure
	double epsilonTP = 0.0;

	/// The second-order partial derivative of the dielectric constant with respect to pressure
	double epsilonPP = 0.0;

	/// The Born function \f$ Z\equiv-\frac{1}{\epsilon} \f$ (see Helgeson and Kirkham, 1974)
	double bornZ = 0.0;

	/// The Born function \f$ Y\equiv\left[\frac{\partial Z}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	double bornY = 0.0;

	/// The Born function \f$ Q\equiv\left[\frac{\partial Z}{\partial P}\right]_{T} \f$ (see Helgeson and Kirkham, 1974)
	double bornQ = 0.0;

	/// The Born function \f$ N\equiv\left[\frac{\partial Q}{\partial P}\right]_{T} \f$ (see Helgeson and Kirkham, 1974)
	double bornN = 0.0;

	/// The Born function \f$ U\equiv\left[\frac{\partial Q}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	double bornU = 0.0;

	/// The Born function \f$ X\equiv\left[\frac{\partial Y}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	double bornX = 0.0;
};

} // namespace Reaktor
