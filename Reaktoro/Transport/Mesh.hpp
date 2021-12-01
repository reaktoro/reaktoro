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

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// A class that defines the mesh for TransportSolver.
class Mesh
{
public:

    /// Construct a Mesh instance.
    Mesh();

    /// Construct a Mesh instance.
    Mesh(Index num_cells, double xl = 0.0, double xr = 1.0);

    /// Initialize a discretization of the mesh.
    auto setDiscretization(Index num_cells, double xl = 0.0, double xr = 1.0) -> void;

    /// Return the number of cells in the mesh.
    auto numCells() const -> Index { return m_num_cells; }

    /// Return the coordinate of the left boundary.
    auto xl() const -> double { return m_xl; }

    /// Return the coordinate of the right boundary (in m).
    auto xr() const -> double { return m_xr; }

    // Return the length of the cells (in m).
    auto dx() const -> double { return m_dx; }

    // Return the x-coordinate of the center of the cells.
    auto xcells() const -> ArrayXrConstRef { return m_xcells; }

private:
    /// The number of cells in the discretization.
    Index m_num_cells = 10;

    /// The x-coordinate of the left boundary (in m).
    double m_xl = 0.0;

    /// The x-coordinate of the right boundary (in m).
    double m_xr = 1.0;

    /// The length of the cells (in m).
    double m_dx = 0.1;

    /// The x-coordinate of the center of the cells.
    ArrayXr m_xcells;
};

} // namespace Reaktoro
