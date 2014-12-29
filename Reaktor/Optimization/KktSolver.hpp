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

struct KktProblem
{
    /// The top-left matrix of the KKT equation
    Matrix H;

    /// The bottom-left matrix of the KKT equation (and its transpose on the top-right corner of the KKT matrix)
    Matrix A;

    /// The top vector of the right-hand side KKT equation
    Vector f;

    /// The bottom vector of the right-hand side KKT equation
    Vector g;
};

struct KktSolution
{
    Vector x;
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

enum KktAlgorithm
{
    KktFullspaceDense     = 0x01,
    KktFullspaceSparse    = 0x02,
    KktRangespace         = 0x04,
    KktNullspace          = 0x08,
    KktNullspacePartial   = 0x10,
};

struct KktOptions
{
    /// The algorithm for the solution of the KKT problem
    KktAlgorithm algorithm = KktFullspaceDense;

    /// The scaling to be applied to the `x` variables
    Vector xscaling;

    /// The scaling to be applied to the `y` variables
    Vector yscaling;

    /// The flag that indicates that the matrix `H` of the KKT equation is diagonal
    bool diagonalH = false;

    /// The flag that indicates that the matrix `A` of the KKT equation is the same for every subsequent KKT calculation
    bool persistentA = false;
};

struct KktResult
{
    KktSolution solution;

    KktStatistics statistics;
};

class KktSolver
{
public:
    KktSolver();

    KktSolver(const KktSolver& other);

    virtual ~KktSolver();

    auto operator=(KktSolver other) -> KktSolver&;

    auto solve(const KktProblem& problem, KktResult& result) -> void;

    auto solve(const KktProblem& problem, KktResult& result, const KktOptions& options) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktor
