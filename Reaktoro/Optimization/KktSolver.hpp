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

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Optimization/Hessian.hpp>

namespace Reaktoro {

/// A type to describe the result of a KKT calculation
struct KktResult
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

/// A type to represent the left-hand side matrix of a KKT equation
/// @see KktSolution, KktVector, KktSolver
struct KktMatrix
{
    /// Construct a custom KktMatrix instance
    KktMatrix(const Hessian& H, const Matrix& A, const Vector& x, const Vector& z)
    : H(H), A(A), x(x), z(z)
    {}

    /// Construct a custom KktMatrix instance
    KktMatrix(const Hessian& H, const Matrix& A, const Vector& x, const Vector& z, double gamma, double delta)
    : H(H), A(A), x(x), z(z), gamma(gamma), delta(delta)
    {}

    /// The Hessian matrix `H` of the KKT matrix equation
    const Hessian& H;

    /// The coefficient matrix `A` of the KKT matrix equation
    const Matrix& A;

    /// The vector of primal variables `x`
    const Vector& x;

    /// The vector of dual variables `z`
    const Vector& z;

    /// The regularization parameter @f$\gamma@f$
    const double gamma = 0.0;

    /// The regularization parameter @f$\delta@f$
    const double delta = 0.0;
};

/// A type to represent the solution vector of a KKT equation
/// @see KktMatrix, KktVector, KktSolver
struct KktSolution
{
    /// The step vector of the primal variables `x`
    Vector dx;

    /// The step vector of the dual variables `y`
    Vector dy;

    /// The step vector of the dual variables `z`
    Vector dz;
};

/// A type to represent the right-hand side vector of a KKT equation
/// @see KktMatrix, KktSolution, KktSolver
struct KktVector
{
    /// The top vector of the right-hand side KKT vector
    Vector rx;

    /// The middle vector of the right-hand side KKT vector
    Vector ry;

    /// The bottom vector of the right-hand side KKT vector
    Vector rz;
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

    /// Return the result of the last calculation
    auto result() const -> const KktResult&;

    /// Set the options for the KKT calculations
    auto setOptions(const KktOptions& options) -> void;

    /// Decompose the KKT matrix before solving it.
    auto decompose(const KktMatrix& lhs) -> void;

    /// Solve the KKT equation using an appropriate and efficient approach
    /// according to a priori decomposition call.
    /// The KKT equation is solved using an LU decomposition if the matrices
    /// `H` and `A` are dense. A rangespace approach is used instead if either
    /// the inverse of the `H` matrix is known or it has a diagonal structure.
    /// Finally, if the matrix `A` has been specified to be constant, then an
    /// efficient nullspace approach is used to reduce the system of linear
    /// equations.
    /// @param rhs The right-hand side vector of the KKT equation
    /// @param sol The solution vector of the KKT equation
    auto solve(const KktVector& rhs, KktSolution& sol) -> void;

private:
    /// Implementation details
    struct Impl;

    /// The pointer to the internal implementation details
    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
