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

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

// Forward declarations
class EquilibriumConstraints;

/// The dimensions of the variables in a constrained equilibrium problem.
/// @see dims(const EquilibriumConstraints&)
struct EquilibriumDims
{
    Index Ne  = 0; ///< The number of chemical elements and electric charge in the chemical system.
    Index Nn  = 0; ///< The number of species in the chemical system.
    Index Npe = 0; ///< The number of equation constraints among the functional constraints.
    Index Npp = 0; ///< The number of property preservation constraints among the functional constraints.
    Index Np  = 0; ///< The number of functional constraints (Np = Npe + Npp).
    Index Nq  = 0; ///< The number of chemical potential constraints.
    Index Nir = 0; ///< The number of reactions prevented from reacting during the equilibrium calculation.
    Index Nc  = 0; ///< The number of introduced control variables (must be equal to Np)
    Index Nx  = 0; ///< The number of variables (Nx = Nn + Np + Nq).
    Index Nb  = 0; ///< The number of components (Nb = Ne + Nir).

    /// Construct a default EquilibriumDims object.
    EquilibriumDims() = default;

    /// Construct an EquilibriumDims object with given equilibrium constraints.
    explicit EquilibriumDims(const EquilibriumConstraints& constraints);
};


} // namespace Reaktoro
