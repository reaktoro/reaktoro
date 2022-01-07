// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "SpeciesThermoProps.hpp"

// Reaktoro includes
#include <Reaktoro/Core/StandardThermoProps.hpp>

namespace Reaktoro {

SpeciesThermoProps::SpeciesThermoProps(const real& T, const real& P, const StandardThermoProps& sprops)
: T(T),
  P(P),
  G0(sprops.G0),
  H0(sprops.H0),
  V0(sprops.V0),
  VT0(sprops.VT0),
  VP0(sprops.VP0),
  Cp0(sprops.Cp0),
  Cv0(VP0 == 0.0 ? Cp0 : Cp0 + T*VT0*VT0/VP0),
  U0(H0 - P*V0),
  S0((H0 - G0)/T),
  A0(U0 - T*S0)
{}

} // namespace Reaktoro
