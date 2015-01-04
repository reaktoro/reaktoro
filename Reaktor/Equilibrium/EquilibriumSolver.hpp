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

// C++ includes
#include <memory>

// Reaktor includes
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

// Forward declarations
class  EquilibriumProblem;
struct EquilibriumOptions;
struct EquilibriumResult;

class EquilibriumSolver
{
public:
    /// Construct a default EquilibriumSolver instance
    EquilibriumSolver();

    /// Construct a copy of an EquilibriumSolver instance
    EquilibriumSolver(const EquilibriumSolver& other);

    /// Destroy this EquilibriumSolver instance
    virtual ~EquilibriumSolver();

    /// Assign a copy of an EquilibriumSolver instance
    auto operator=(EquilibriumSolver other) -> EquilibriumSolver&;

    /// Find a initial guess for an equilibrium problem
    /// @param problem The definition of the equilibrium problem
    /// @param result[in,out] The initial guess and the final result of the equilibrium approximation
    /// @param options The options for the equilibrium calculation
    auto approximate(const EquilibriumProblem& problem, EquilibriumResult& result) -> void;

    /// Find a initial guess for an equilibrium problem with given options
    /// @param problem The definition of the equilibrium problem
    /// @param result[in,out] The initial guess and the final result of the equilibrium approximation
    /// @param options The options for the equilibrium calculation
    auto approximate(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void;

    /// Solve an equilibrium problem
    /// @param problem The definition of the equilibrium problem
    /// @param result[in,out] The initial guess and the final result of the equilibrium calculation
    /// @param options The options for the equilibrium calculation
    auto solve(const EquilibriumProblem& problem, EquilibriumResult& result) -> void;

    /// Solve an equilibrium problem with given options
    /// @param problem The definition of the equilibrium problem
    /// @param result[in,out] The initial guess and the final result of the equilibrium calculation
    /// @param options The options for the equilibrium calculation
    auto solve(const EquilibriumProblem& problem, EquilibriumResult& result, const EquilibriumOptions& options) -> void;

    /// Compute the partial derivatives @f$\left.\frac{\partial n}{\partial T}\right|_{P,b}@f$.
    /// These derivatives tell us how much the equilibrium composition
    /// of the equilibrium species will change with an infinitesimal change
    /// in temperature. They are useful when solving non-linear problems that
    /// involve equilibrium calculations and derivatives with respect to temperature.
    /// @param result The result of an equilibrium calculation performed a priori
    auto dndt(const EquilibriumResult& result) -> Vector;

    /// Compute the partial derivatives @f$\left.\frac{\partial n}{\partial P}\right|_{T,b}@f$.
    /// These derivatives tell us how much the equilibrium composition
    /// of the equilibrium species will change with an infinitesimal change
    /// in pressure. They are useful when solving non-linear problems that
    /// involve equilibrium calculations and derivatives with respect to pressure.
    /// @param result The result of an equilibrium calculation performed a priori
    auto dndp(const EquilibriumResult& result) -> Vector;

    /// Compute the partial derivatives @f$\left.\frac{\partial n}{\partial b}\right|_{T,P}@f$.
    /// These derivatives tell us how much the equilibrium composition
    /// of the equilibrium species will change with an infinitesimal change
    /// in the amounts of elements. They are useful when solving non-linear problems that
    /// involve equilibrium calculations and derivatives with respect to element amounts.
    /// @param result The result of an equilibrium calculation performed a priori
    auto dndb(const EquilibriumResult& result) -> Matrix;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
