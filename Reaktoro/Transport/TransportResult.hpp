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

namespace Reaktoro {

/// Provide timing information of the operations in a transport time step calculation.
struct TransportTiming
{
    /// The time spent during the last time step of a transport calculation.
    double step = 0.0;

    /// The time spent during the last time step for matrix equation assembly.
    double matrix_equation_assembly = 0.0;

    /// The time spent during the last time step for solving the discrete transport matrix equation.
    double matrix_equation_solve = 0.0;

    /// Self addition of another TransportTiming instance to this one.
    auto operator+=(const TransportTiming& other) -> TransportTiming&;
};

/// Provide result information of a transport time step calculation.
struct TransportResult
{
    /// The timing information of the operations in a transport time step calculation.
    TransportTiming timing;
};

} // namespace Reaktoro
