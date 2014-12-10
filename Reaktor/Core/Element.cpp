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

#include "Element.hpp"

namespace Reaktor {

struct Element::Impl
{
    ElementData data;
};

Element::Element()
: pimpl(new Impl())
{}

Element::Element(const ElementData& data)
: pimpl(new Impl())
{
    pimpl->data = data;
}

auto Element::name() const -> std::string
{
    return pimpl->data.name;
}

auto Element::molarMass() const -> double
{
    return pimpl->data.molar_mass;
}

auto operator<(const Element& lhs, const Element& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

} // namespace Reaktor
