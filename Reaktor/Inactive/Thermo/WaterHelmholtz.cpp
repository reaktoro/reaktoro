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

#include "WaterHelmholtz.hpp"

// Reaktor includes
#include <Reaktor/Thermo/WaterThermoWagnerPruss.hpp>
#include <Reaktor/Thermo/WaterThermoHGK.hpp>

namespace Reaktor {

WaterHelmholtz::WaterHelmholtz()
: helmholtz(0), helmholtzT(0), helmholtzD(0), helmholtzTT(0), helmholtzTD(0),
  helmholtzDD(0), helmholtzTTT(0), helmholtzTTD(0), helmholtzTDD(0), helmholtzDDD(0)
{}

auto waterHelmholtzState(double T, double D) -> WaterHelmholtz
{
	return waterHelmholtzState(T, D, WagnerPruss);
}

auto waterHelmholtzState(double T, double D, WaterThermoModel model) -> WaterHelmholtz
{
	switch(model)
	{
	case WagnerPruss: return waterHelmholtzStateWagnerPruss(T, D);
	case HGK:         return waterHelmholtzStateHGK(T, D);
	default:          return waterHelmholtzStateWagnerPruss(T, D);
	}
}

} // namespace Reaktor
