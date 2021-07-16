// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <vector>

namespace Reaktoro {

template<typename T> using Table1D = std::vector<T>;
template<typename T> using Table2D = std::vector<std::vector<T>>;
template<typename T> using Table3D = std::vector<std::vector<std::vector<T>>>;

template<typename T>
auto table1D(unsigned dim1) -> Table1D<T>
{
    return Table1D<T>(dim1);
}

template<typename T>
auto table2D(unsigned dim1, unsigned dim2) -> Table2D<T>
{
    return Table2D<T>(dim1, table1D<T>(dim2));
}

template<typename T>
auto table3D(unsigned dim1, unsigned dim2, unsigned dim3) -> Table3D<T>
{
    return Table3D<T>(dim1, table2D<T>(dim2, dim3));
}

} // namespace Reaktoro
