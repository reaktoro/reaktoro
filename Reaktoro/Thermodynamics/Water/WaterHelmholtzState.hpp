// Reaktoro is a C++ library for computational reaction modelling.
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
