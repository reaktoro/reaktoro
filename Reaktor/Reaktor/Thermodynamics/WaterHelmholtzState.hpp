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

struct WaterHelmholtzState
{
	/// The specific Helmholtz free energy of water (in units of J/kg)
	double helmholtz = 0.0;

	/// The first-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	double helmholtzT = 0.0;

	/// The first-order partial derivative of the specific Helmholtz free energy of water with respect to density
	double helmholtzD = 0.0;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	double helmholtzTT = 0.0;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to temperature and density
	double helmholtzTD = 0.0;

	/// The second-order partial derivative of the specific Helmholtz free energy of water with respect to density
	double helmholtzDD = 0.0;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature
	double helmholtzTTT = 0.0;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature, temperature, and density
	double helmholtzTTD = 0.0;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to temperature, density, and density
	double helmholtzTDD = 0.0;

	/// The third-order partial derivative of the specific Helmholtz free energy of water with respect to density
	double helmholtzDDD = 0.0;
};

} // namespace Reaktor
