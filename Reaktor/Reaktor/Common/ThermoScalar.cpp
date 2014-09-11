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

#include "ThermoScalar.hpp"

// Reaktor includes
#include <Reaktor/Common/ThermoVector.hpp>

namespace Reaktor {

ThermoScalar::ThermoScalar()
: val(), ddt(), ddp(), ddn()
{}

ThermoScalar::ThermoScalar(double val, double ddt, double ddp, const Vector& ddn)
: val(val), ddt(ddt), ddp(ddp), ddn(ddn)
{}

ThermoScalar::ThermoScalar(unsigned ndim)
: val(), ddt(), ddp(), ddn(ndim)
{}

ThermoScalar::ThermoScalar(const ThermoVectorRow& row)
{
    val = row.val[0];
    ddt = row.ddt[0];
    ddp = row.ddp[0];
    ddn = row.ddn;
}

auto ThermoScalar::zero(unsigned nspecies) -> ThermoScalar
{
    ThermoScalar scalar;
    scalar.val = 0.0;
    scalar.ddt = 0.0;
    scalar.ddp = 0.0;
    scalar.ddn = zeros(nspecies);
    return scalar;
}

} // namespace Reaktor
