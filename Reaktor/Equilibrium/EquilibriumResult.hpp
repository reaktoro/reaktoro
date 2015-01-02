// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Optimization/OptimumResult.hpp>

namespace Reaktor {

/// A type used to describe the solution of an equilibrium calculation
/// @see EquilibriumStatistics, EquilibriumResult
struct EquilibriumSolution
{
    /// The molar amounts of the species (in units of mol)
    Vector n;

	/// The chemical potentials of the species (in units of J/mol)
	Vector u;

    /// The molar amounts of the equilibrium species (in units of mol)
    Vector ne;

	/// The Lagrange multipliers w.r.t. the equilibrium elements and charge (in units of J/mol)
	Vector ye;

	/// The Lagrange multipliers w.r.t. the equilibrium species (in units of J/mol)
	Vector ze;

	/// The Gibbs energy of the equilibrium partition at equilibrium (in units of J)
	double ge;

	/// The chemical potentials of the equilibrium species (in units of J/mol)
	Vector ue;

	/// The partial derivative of the equilibrium species abundances w.r.t. temperature (in units of mol/K)
	Vector dndt;

	/// The partial derivative of the equilibrium species abundances w.r.t. pressure (in units of mol/Pa)
	Vector dndp;

	/// The partial derivative of the equilibrium species abundances w.r.t. equilibrium element abundances (in units of mol/mol)
	Matrix dndb;

	/// The indices of the stable equilibrium species
	Indices indices_stable_species;
};

/// A type used to describe the statistics of an equilibrium calculation
/// @see EquilibriumSolution, EquilibriumResult
struct EquilibriumStatistics : OptimumStatistics
{
    /// Construct a default EquilibriumStatistics instance
    EquilibriumStatistics();

    /// Construct an EquilibriumStatistics instance from a OptimumStatistics instance
    EquilibriumStatistics(const OptimumStatistics& other);
};

/// A type used to describe the result of an equilibrium calculation, with its solution and statistics
/// @see EquilibriumSolution. EquilibriumStatistics
struct EquilibriumResult
{
    /// The solution of the equilibrium calculation
	EquilibriumSolution solution;

    /// The statistics of the equilibrium calculation
	EquilibriumStatistics statistics;

	/// The result of the optimisation calculation
	OptimumResult optimum;
};

} // namespace Reaktor
