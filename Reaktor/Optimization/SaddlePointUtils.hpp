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
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

struct SaddlePointProblem
{
    Matrix H;
    Matrix A;
    Vector f;
    Vector g;
};

struct SaddlePointSolution
{
    Vector x;
    Vector y;
};

struct SaddlePointStatistics
{
    /// The flag that indicates if the saddle point calculation converged
    bool converged = false;

    /// The number of iterations for the solution of the saddle point problem
    unsigned num_iterations = 0;

    /// The wall time spent for the solution of the saddle point problem (in units of s)
    double time = 0;
};

struct SaddlePointInternal
{
    Matrix Z;
    Matrix Y;
    Matrix ZtHZ;
    Vector xZ;
    Matrix L;
    Matrix U;
    Matrix P;
    Matrix R;
    Matrix lhs;
    Vector rhs;
    Vector sol;
};

enum SaddlePointAlgorithm
{
    FullspaceDense     = 0x01,
    FullspaceSparse    = 0x02,
    Rangespace         = 0x04,
    Nullspace          = 0x08,
    NullspacePartial   = 0x10,
};

enum SaddlePointProperties
{
    DiagonalH = 0x01,
    ConstantA = 0x02
};

struct SaddlePointOptions
{
    SaddlePointAlgorithm algorithm = FullspaceDense;

    SaddlePointProperties properties;
};

struct SaddlePointResult
{
    SaddlePointSolution solution;

    SaddlePointStatistics statistics;

    SaddlePointInternal internal;
};

auto solve(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void;

auto solveFullspaceDense(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void;

auto solveFullspaceSparse(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void;

auto solveRangespace(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void;

auto solveNullspace(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void;

auto solveNullspacePartial(const SaddlePointProblem& problem, SaddlePointResult& result, const SaddlePointOptions& options) -> void;

/// Calculate the null and range space of a matrix.
/// This method calculates the null and range space matrices **Z** and **Y**
/// such that **AZ = 0** and **AY = I**, where **I** is the identity matrix.
/// @param A The matrix whose range and null space is to be calculated
/// @param[out] Z The null space matrix
/// @param[out] Y The range space matrix
auto nullrange(const Matrix& A, Matrix& Z, Matrix& Y) -> void;

} // namespace Reaktor
