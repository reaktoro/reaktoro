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

// C++ includes
#include <memory>

namespace Reaktoro {

//// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumState;

/// A type that represents a regularized optimization problem.
class Regularizer
{
public:
    /// Construct a default Regularizer instance.
    Regularizer();

    /// Construct a copy of an Regularizer instance
    Regularizer(const Regularizer& other);

    /// Destroy this instance
    virtual ~Regularizer();

    /// Assign an Regularizer instance to this instance
    auto operator=(Regularizer other) -> Regularizer&;

    /// Regularize the optimum problem, state, and options before they are used in an optimization calculation.
    /// @param problem The optimum problem to be regularized.
    /// @param state The optimum state to be regularized.
    /// @param options The optimum options to be regularized.
    auto regularize(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void;

    /// Recover an optimum state to an state that corresponds to the original optimum problem.
    /// @param problem The optimum problem regularized in method `regularize`.
    /// @param state[in,out] The optimum state regularized in method `regularize`.
    auto recover(const OptimumProblem& problem, OptimumState& state) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
