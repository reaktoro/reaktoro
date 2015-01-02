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

/// A type to describe the statistics of the solution of a KKT equation
struct KktStatistics
{
    /// The flag that indicates if the KKT calculation succeeded
    bool succeeded = false;

    /// The wall time spent for the solution of the KKT problem (in units of s)
    double time = 0;
};

/// A enumeration of possible schemes for the solution of the KKT equation
enum class KktScheme
{
    Full, InverseH, DiagonalH, ConstantA, Uninitialised
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

    /// Return the statistics of the last KKT calculation.
    auto statistics() const -> const KktStatistics&;

    /// Return the active scheme for the solution of a KKT equation
    auto scheme() const -> KktScheme;

    /// Specify to the solver that the matrix `A` is constant.
    /// Calling this method ensures that the nullspace and rangespace
    /// matrices `Z` and `Y` (i.e., `AZ = 0` and `AY = I`) are computed
    /// only once and reused later when solving the KKT equation.
    /// @see decomposeWithConstantA
    auto setConstantA(const Matrix& A) -> void;

    /// Pre-decompose the KKT matrix using dense matrices.
    /// @param H The dense top-left corner matrix of the KKT equation
    /// @param A The dense bottom-left corner matrix of the KKT equation
    auto decompose(const Matrix& H, const Matrix& A) -> void;

    /// Pre-decompose the KKT matrix taking advantage of known inverse matrix `inv(H)`.
    /// @param invH The inverse of the dense top-left corner matrix of the KKT equation
    /// @param A The dense bottom-left corner matrix of the KKT equation
    auto decomposeWithInverseH(const Matrix& invH, const Matrix& A) -> void;

    /// Pre-decompose the KKT matrix taking advantage of the diagonal structure of the matrix `H`.
    /// @param H The diagonal top-left corner matrix (represented as a vector here) of the KKT equation
    /// @param A The dense bottom-left corner matrix of the KKT equation
    auto decomposeWithDiagonalH(const Vector& H, const Matrix& A) -> void;

    /// Pre-decompose the KKT matrix taking advantage of the unchanged behaviour of the matrix `A`.
    /// Note that the method `setConstantA` must be called before to set the `A` matrix.
    /// @param H The dense top-left corner matrix of the KKT equation
    /// @see setConstantA
    auto decomposeWithConstantA(const Matrix& H) -> void;

    /// Solve the KKT equation using an appropriate and efficient approach according to a priori decomposition call.
    /// The KKT equation is solved using an LU decomposition if the matrices `H` and `A` are dense.
    /// A rangespace approach is used instead if either the inverse of the `H` matrix is known or it
    /// has a diagonal structure. Finally, if the matrix `A` has been specified to be constant,
    /// then an efficient nullspace approach is used to reduce the system of linear equations.
    /// @param a The top vector of the right-hand side of the KKT equation
    /// @param b The bottom vector of the right-hand side of the KKT equation
    /// @param x The top vector of unknowns of the KKT equation
    /// @param y The bottom vector of unknowns of the KKT equation
    auto solve(const Vector& a, const Vector& b, Vector& x, Vector& y) -> void;

private:
    /// Implementation details
    struct Impl;

    /// The pointer to the internal implementation details
    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
