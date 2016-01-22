// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// Reaktoro includes
#include<Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

struct WaterElectroState
{
	/// The dielectric constant of water
	ThermoScalar epsilon;

	/// The first-order partial derivative of the dielectric constant with respect to temperature
	ThermoScalar epsilonT;

	/// The first-order partial derivative of the dielectric constant with respect to pressure
	ThermoScalar epsilonP;

	/// The second-order partial derivative of the dielectric constant with respect to temperature
	ThermoScalar epsilonTT;

	/// The second-order partial derivative of the dielectric constant with respect to temperature and pressure
	ThermoScalar epsilonTP;

	/// The second-order partial derivative of the dielectric constant with respect to pressure
	ThermoScalar epsilonPP;

	/// The Born function \f$ Z\equiv-\frac{1}{\epsilon} \f$ (see Helgeson and Kirkham, 1974)
	ThermoScalar bornZ;

	/// The Born function \f$ Y\equiv\left[\frac{\partial Z}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	ThermoScalar bornY;

	/// The Born function \f$ Q\equiv\left[\frac{\partial Z}{\partial P}\right]_{T} \f$ (see Helgeson and Kirkham, 1974)
	ThermoScalar bornQ;

	/// The Born function \f$ N\equiv\left[\frac{\partial Q}{\partial P}\right]_{T} \f$ (see Helgeson and Kirkham, 1974)
	ThermoScalar bornN;

	/// The Born function \f$ U\equiv\left[\frac{\partial Q}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	ThermoScalar bornU;

	/// The Born function \f$ X\equiv\left[\frac{\partial Y}{\partial T}\right]_{P} \f$ (see Helgeson and Kirkham, 1974)
	ThermoScalar bornX;
};

} // namespace Reaktoro
