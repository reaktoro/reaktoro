// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "Class.hpp"

namespace Reaktoro {

struct Class::Impl
{

};

Class::Class()
: pimpl(new Impl())
{
}

Class::Class(const Class& other)
: pimpl(new Impl(*other.pimpl))
{
}

Class::~Class()
{
}

auto Class::operator=(Class other) -> Class&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

} // namespace Reaktoro
