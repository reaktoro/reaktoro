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

struct KktSolution
{
    /// The `x` variables of the KKT equation
    Vector x;

    /// The `y` variables of the KKT equation
    Vector y;
};

struct KktStatistics
{
    /// The flag that indicates if the saddle point calculation converged
    bool converged = false;

    /// The number of iterations for the solution of the saddle point problem
    unsigned num_iterations = 0;

    /// The wall time spent for the solution of the saddle point problem (in units of s)
    double time = 0;
};

struct KktResult
{
    KktSolution solution;

    KktStatistics statistics;
};

struct KktMatrix
{
    /// The top-left matrix `H` of the KKT equation
    Matrix H;

    /// The bottom-left matrix `A` of the KKT equation
    Matrix A;

    /// The inverse of the top-left matrix `H` of the KKT equation.
    /// If this inverse is provided, then an appropriate algorithm based
    /// on a rangespace approach will be used to solve the KKT equation.
    Matrix invH;

    /// The top-left matrix `H` of the KKT equation as a diagonal matrix
    /// If the `H` matrix is diagonal, then an appropriate algorithm based
    /// on a rangespace approach will be used to solve the KKT equation.
    Vector diagH;
};

struct KktVector
{
    /// The top vector of the right-hand side KKT equation
    Vector f;

    /// The bottom vector of the right-hand side KKT equation
    Vector g;
};

struct KktScaling
{
    /// The scaling to be applied to the `x` variables
    Vector x;

    /// The scaling to be applied to the `y` variables
    Vector y;
};

class KktProblem
{
public:
    KktProblem();

    KktProblem(const KktProblem& other);

    virtual ~KktProblem();

    auto operator=(KktProblem other) -> KktProblem&;

    auto setConstantA(const Matrix& A) -> void;

    auto decompose(const KktMatrix& lhs) -> void;

    auto decomposeWithInverseH(const KktMatrix& lhs) -> void;

    auto decomposeWithDiagonalH(const KktMatrix& lhs) -> void;

    auto decomposeWithConstantA(const KktMatrix& lhs) -> void;

    auto solve(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void;

    auto solveWithInverseH(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void;

    auto solveWithDiagonalH(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void;

    auto solveWithConstantA(const KktMatrix& lhs, const KktVector& rhs, KktResult& result) -> void;

private:
    /// Implementation details
    struct Impl;

    /// The pointer to the internal implementation details
    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
