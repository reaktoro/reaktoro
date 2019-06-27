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

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// A class that defines the mesh for TransportSolver.
class Mesh
{
public:
    Mesh();

    Mesh(Index num_cells, double xl = 0.0, double xr = 1.0);

    auto setDiscretization(Index num_cells, double xl = 0.0, double xr = 1.0) -> void;

    auto numCells() const -> Index { return m_num_cells; }

    auto xl() const -> double { return m_xl; }

    auto xr() const -> double { return m_xr; }

    auto dx() const -> double { return m_dx; }

    auto xcells() const -> VectorConstRef { return m_xcells; }

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
    Vector m_xcells;
};

} // namespace Reaktoro
