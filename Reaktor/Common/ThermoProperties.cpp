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

#include "ThermoProperties.hpp"

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>

namespace Reaktor {
namespace internal {

auto assertThermoProperties(const Vector& val, const Vector& ddt, const Vector& ddp) -> void
{
    Assert(val.size() == ddt.size() and val.size() == ddt.size(),
        "ThermoProperties requires arguments with the same dimensions.");
}

} // namespace internal

ThermoProperties::ThermoProperties(const Vector& val, const Vector& ddt, const Vector& ddp)
: m_val(val), m_ddt(ddt), m_ddp(ddp)
{
    internal::assertThermoProperties(val, ddt, ddp);
}

ThermoProperties::ThermoProperties(Vector&& val, Vector&& ddt, Vector&& ddp)
{
    internal::assertThermoProperties(val, ddt, ddp);
}

auto ThermoProperties::val() const -> Vector
{
	return m_val;	
}

auto ThermoProperties::ddt() const -> Vector
{
	return m_ddt;	
}

auto ThermoProperties::ddp() const -> Vector
{
	return m_ddp;	
}

}  // namespace Reaktor
