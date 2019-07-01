// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <ostream>
#include <vector>

namespace Reaktoro {

/// Output the contents of a std::vector container to the output stream.
template<typename T>
auto operator<<(std::ostream& out, const std::vector<T>& items) -> std::ostream&
{
    for(std::size_t i = 0; i < items.size(); ++i)
        out << (i > 0 ? ", " : "") << item[i];
    out << std::endl;
}

} // namespace Reaktoro
