// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
class EquilibriumSpecs;

/// The dimensions of the variables in a chemical equilibrium problem.
///
/// In a chemical equilibrium problem, we are interested in computing the
/// values of the following variables:
///
/// * @eq{n}, the vector with the amounts of the chemical species;
/// * @eq{p}, the vector with the *p* control variables;
/// * @eq{q}, the vector with the *q* control variables;
///
/// The *p* control variables are those variables introduced to enforce
/// constraints of *equation type*. These are *temperature* and/or *pressure*, in
/// case they are unknown, as well as the amounts of *explicit titrants* (the
/// substances that we permit in/out of the system).
///
/// The *q* control variables are those variables introduced to enforce
/// constraints of *chemical potential type*. These are *implicit titrants*
/// (also substances allowed in/out of the system) whose amounts need to be
/// found so that the chemical potential constraints are satisfied.
struct EquilibriumDims
{
    Index Ne = 0; ///< The number of elements in the chemical system.
    Index Nb = 0; ///< The number of components in the chemical equilibrium problem (`Nb = Ne` if no inert reactions are specified).
    Index Nn = 0; ///< The number of species in the chemical system.
    Index Np = 0; ///< The number of *p* control variables (temperature, pressure, and amounts of explicit titrants when these are introduced unknowns).
    Index Nq = 0; ///< The number of *q* control variables (the amounts of the implicit titrants when these are introduced unknowns).
    Index Nt = 0; ///< The number of substances for which the chemical system is open to (the number of explicit and implicit titrants).
    Index Nx = 0; ///< The number of variables *x* in *x = (n, q)* (equivalent to `Nn + Nq`).
    Index Nu = 0; ///< The number of unknown variables in the chemical equilibrium problem (equivalent to `Nn + Np + Nq`).

    /// Construct a default EquilibriumDims object.
    EquilibriumDims() = default;

    /// Construct an EquilibriumDims object with given equilibrium specifications.
    explicit EquilibriumDims(const EquilibriumSpecs& specs);
};

} // namespace Reaktoro
