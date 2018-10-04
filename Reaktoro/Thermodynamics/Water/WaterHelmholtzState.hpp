// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

struct WaterHelmholtzState
{
	/// The specific Helmholtz free energy of water (in units of J/kg)
	ThermoScalar helmholtz;

	/// The first-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	ThermoScalar helmholtzT;

	/// The first-order partial derivative of the specific Helmholtz free energy of water with respect to density
	ThermoScalar helmholtzD;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	ThermoScalar helmholtzTT;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to temperature and density
	ThermoScalar helmholtzTD;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to density
	ThermoScalar helmholtzDD;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	ThermoScalar helmholtzTTT;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature, temperature, and density
	ThermoScalar helmholtzTTD;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature, density, and density
	ThermoScalar helmholtzTDD;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to density
	ThermoScalar helmholtzDDD;
};

} // namespace Reaktoro
