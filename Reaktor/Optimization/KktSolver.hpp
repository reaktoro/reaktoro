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
struct OptimumState;

/// A type to describe some information about the KKT calculation
struct KktInfo
{
    /// The flag that indicates if the KKT calculation succeeded
    bool succeeded = false;

    /// The wall time spent for the decomposition of the KKT problem (in units of s)
    double time_decompose = 0;

    /// The wall time spent for the solution of the KKT problem (in units of s)
    double time_solve = 0;
};

/// An enumeration of possible methods for the solution of a KKT equation
enum class KktMethod
{
    /// Use a partial pivoting LU algorithm on the full KKT equation.
    /// This can only be used for dense Hessian matrices.
    PartialPivLU,

    /// Use a full pivoting LU algorithm on the full KKT equation.
    /// This can only be used for dense Hessian matrices.
    FullPivLU,

    /// Use a nullspace method to solve the KKT equation.
    /// This method is advisable when there are many equality constraints
    /// and these are linear so that their gradient is constant.
    Nullspace,

    /// Use a rangespace method to solve the KKT equation.
    /// This method is advisable when the Hessian matrix can be easily
    /// inverted such as a quasi-Newton approximation or a diagonal matrix.
    Rangespace,

    /// Use a method that fits better to the type of KKT equation.
    /// This option will ensure that a rangespace method is used when
    /// the Hessian matrix is diagonal or its inverse is available.
    /// It will use a `PartialPivLU` method for dense KKT equations.
    Automatic,
};

/// A type to describe the options for the KKT calculation
struct KktOptions
{
    /// The method for the solution of the KKT equations
    KktMethod method = KktMethod::Automatic;
};

/// A type to describe a solver for a KKT equation
class KktSolver
{
public:
    /// Construct a default KktSolver instance
    KktSolver();

    /// Construct a copy of a KktSolver instance
    KktSolver(const KktSolver& other);

    /// Destroy this KktSolver instance
    virtual ~KktSolver();

    /// Assign a KktSolver instance to this
    auto operator=(KktSolver other) -> KktSolver&;

    /// Return the info of the last calculation
    auto info() const -> const KktInfo&;

    /// Set the options for the KKT calculations
    auto setOptions(const KktOptions& options) -> void;

    /// Decompose the KKT matrix before solving it.
    auto decompose(const OptimumState& state) -> void;

    /// Solve the KKT equation using an appropriate and efficient approach
    /// according to a priori decomposition call.
    /// The KKT equation is solved using an LU decomposition if the matrices
    /// `H` and `A` are dense. A rangespace approach is used instead if either
    /// the inverse of the `H` matrix is known or it has a diagonal structure.
    /// Finally, if the matrix `A` has been specified to be constant, then an
    /// efficient nullspace approach is used to reduce the system of linear
    /// equations.
    /// @param a The top vector of the right-hand side of the KKT equation
    /// @param b The bottom vector of the right-hand side of the KKT equation
    /// @param dx The step on the primal variables `x`
    /// @param dy The step on the dual variables `y`
    auto solve(const Vector& a, const Vector& b, Vector& dx, Vector& dy) -> void;

private:
    /// Implementation details
    struct Impl;

    /// The pointer to the internal implementation details
    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
