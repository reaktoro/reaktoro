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
};

struct SaddlePointResult
{
    SaddlePointSolution solution;
    SaddlePointInternal internal;
};

auto solveFull(const SaddlePointProblem& problem, SaddlePointResult& result) -> void;

auto solveNullspace(const SaddlePointProblem& problem, SaddlePointResult& result) -> void;

} // namespace Reaktor
