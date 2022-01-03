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

#include "NasaThermoParams.hpp"

namespace Reaktoro {

auto operator!=(const NasaThermoParams& l, const NasaThermoParams& r) -> bool
{
    return
        l.Tmin != r.Tmin &&
        l.Tmax != r.Tmax &&
        l.qN != r.qN &&
        l.q1 != r.q1 &&
        l.q2 != r.q2 &&
        l.q3 != r.q3 &&
        l.q4 != r.q4 &&
        l.q5 != r.q5 &&
        l.q6 != r.q6 &&
        l.q7 != r.q7 &&
        l.a1 != r.a1 &&
        l.a2 != r.a2 &&
        l.a3 != r.a3 &&
        l.a4 != r.a4 &&
        l.a5 != r.a5 &&
        l.a6 != r.a6 &&
        l.a7 != r.a7 &&
        l.b1 != r.b1 &&
        l.b2 != r.b2;
}

auto operator==(const NasaThermoParams& l, const NasaThermoParams& r) -> bool
{
    return !(l != r);
}

} // namespace Reaktoro
