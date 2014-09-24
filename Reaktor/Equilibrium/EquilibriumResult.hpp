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
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

/// A type used to describe the solution of an equilibrium calculation
/// @see EquilibriumStatistics, EquilibriumResult
struct EquilibriumSolution
{
    /// The abundance of the chemical species (in units of mol)
	Vector n;

	/// The Lagrange multipliers w.r.t. the chemical elements and charge (in units of J/mol)
	Vector y;

	/// The Lagrange multipliers w.r.t. the chemical species (in units of J/mol)
	Vector z;

	/// The partial derivative of the species abundances w.r.t. temperature (in units of mol/K)
	Vector dndt;

	/// The partial derivative of the species abundances w.r.t. pressure (in units of mol/Pa)
	Vector dndp;

	/// The partial derivative of the species abundances w.r.t. element abundances (in units of mol/mol)
	Matrix dndb;

	/// The indices of the stable species at equilibrium
	Indices indices_stable_species;
};

/// A type used to describe the statistics of an equilibrium calculation
/// @see EquilibriumSolution, EquilibriumResult
struct EquilibriumStatistics
{
    /// The flag that indicates if the equilibrium calculation converged
	bool converged;

	/// The number of iterations executed in the equilibrium calculation
	unsigned num_iterations;

	/// The number of function evaluations (i.e., Gibbs energy function) executed in the equilibrium calculation
	unsigned num_func_evals;

	/// The number of gradient evaluations executed in the equilibrium calculation
	unsigned num_grad_evals;

	/// The number of Hessian evaluations executed in the equilibrium calculation
	unsigned num_hessian_evals;
};

/// A type used to describe the internal data used in an equilibrium calculation
struct EquilibriumInternal
{

};

/// A type used to describe the result of an equilibrium calculation, with its solution and statistics
/// @see EquilibriumSolution. EquilibriumStatistics
struct EquilibriumResult
{
    /// The solution of the equilibrium calculation
	EquilibriumSolution solution;

    /// The statistics of the equilibrium calculation
	EquilibriumStatistics statistics;

    /// The internal data of the equilibrium calculation
	EquilibriumInternal internal;
};

} // namespace Reaktor
