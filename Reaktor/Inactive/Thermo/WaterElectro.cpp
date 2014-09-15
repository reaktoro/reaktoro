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

#include "WaterElectro.hpp"

// C++ includes
#include <cmath>
using namespace std;

// Reaktor includes
#include <Reaktor/Thermo/WaterElectroJohnsonNorton.hpp>
#include <Reaktor/Thermo/WaterThermo.hpp>

namespace Reaktor {

WaterElectro::WaterElectro()
: epsilon(0), epsilonT(0), epsilonP(0), epsilonTT(0), epsilonTP(0), epsilonPP(0),
  bornZ(0), bornY(0), bornQ(0), bornN(0), bornU(0), bornX(0)
{}

auto operator<<(std::ostream& out, const WaterElectro& we) -> std::ostream&
{
	out << "epsilon   = " << we.epsilon   << std::endl;
	out << "epsilonT  = " << we.epsilonT  << std::endl;
	out << "epsilonP  = " << we.epsilonP  << std::endl;
	out << "epsilonTT = " << we.epsilonTT << std::endl;
	out << "epsilonTP = " << we.epsilonTP << std::endl;
	out << "epsilonPP = " << we.epsilonPP << std::endl;
	out << "bornZ     = " << we.bornZ     << std::endl;
	out << "bornY     = " << we.bornY     << std::endl;
	out << "bornQ     = " << we.bornQ     << std::endl;
	out << "bornN     = " << we.bornN     << std::endl;
	out << "bornU     = " << we.bornU     << std::endl;
	out << "bornX     = " << we.bornX     << std::endl;

	return out;
}

auto waterElectroState(double T, double P, const WaterThermo& wt) -> WaterElectro
{
	return waterElectroStateJohnsonNorton(T, P, wt);
}

} // namespace Reaktor
