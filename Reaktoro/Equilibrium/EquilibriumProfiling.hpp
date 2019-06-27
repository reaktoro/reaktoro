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

namespace Reaktoro {

/// Provide profiling information of the operations during an equilibrium calculation.
struct EquilibriumProfiling
{
    /// The time spent for solving the chemical equilibrium problem.
    double time_solve = 0.0;

    /// The time spent for computing the standard thermochemical properties of the system.
    double time_standard_properties = 0.0;

    /// The time spent for computing the chemical properties of the system.
    double time_chemical_properties = 0.0;

    /// The time spent for computing the sensitivity derivatives of the chemical state.
    double time_sensitivity = 0.0;
};

} // namespace Reaktoro
