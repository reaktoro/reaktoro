// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Optimization/OptimumSolverBase.hpp>

namespace Reaktoro {

// Forward declarations
struct OptimumProblem;
struct OptimumResult;
struct OptimumState;
struct OptimumOptions;

class OptimumSolverIpAction : public OptimumSolverBase
{
public:
    OptimumSolverIpAction();

    OptimumSolverIpAction(const OptimumSolverIpAction& other);

    virtual ~OptimumSolverIpAction();

    auto operator=(OptimumSolverIpAction other) -> OptimumSolverIpAction&;

    virtual auto solve(const OptimumProblem& problem, OptimumState& state) -> OptimumResult;

    virtual auto solve(const OptimumProblem& problem, OptimumState& state, const OptimumOptions& options) -> OptimumResult;

    virtual auto clone() const -> OptimumSolverBase*;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
