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

#pragma once

// C++ includes
#include <iostream>

namespace Reaktoro {

/// The list of states of matter.
enum class StateOfMatter
{
    Solid, Liquid, Gas, Plasma
};

/// Output a StateOfMatter value.
inline auto operator<<(std::ostream& out, StateOfMatter option) -> std::ostream&
{
    switch(option)
    {
    case StateOfMatter::Solid:  out << "Solid";  return out;
    case StateOfMatter::Liquid: out << "Liquid"; return out;
    case StateOfMatter::Gas:    out << "Gas";    return out;
    case StateOfMatter::Plasma: out << "Plasma"; return out;
    }
}

} // namespace Reaktoro
