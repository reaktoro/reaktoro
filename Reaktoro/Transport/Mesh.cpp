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

#include "Mesh.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

Mesh::Mesh()
{}

Mesh::Mesh(Index num_cells, double xl, double xr)
{
    setDiscretization(num_cells, xl, xr);
}

auto Mesh::setDiscretization(Index num_cells, double xl, double xr) -> void
{
    Assert(xr > xl, "Could not set the discretization.",
        "The x-coordinate of the right boundary needs to be "
        "larger than that of the left boundary.");

    m_num_cells = num_cells;
    m_xl = xl;
    m_xr = xr;
    m_dx = (xr - xl) / num_cells;
    m_xcells = linspace(num_cells, xl + 0.5*m_dx, xr - 0.5*m_dx);
}

} // namespace Reaktoro
