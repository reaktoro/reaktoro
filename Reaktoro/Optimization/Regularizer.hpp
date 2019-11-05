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

namespace Reaktoro {

//// Forward declarations
struct OptimumOptions;
struct OptimumProblem;
struct OptimumState;

/// A type that describes the options for regularizing linear constraints.
struct RegularizerOptions
{
    /// The boolean flag that indicates if echelonization should be performed.
    /// The echelonization of the linear constraints can help on robustness and
    /// accuracy by minimizing round-off errors.
    bool echelonize = true;

    /// The maximum denominator that can exist in the coefficient matrix `A`.
    /// Set this option to zero if the coefficients in `A` are not represented
    /// by rational numbers. Otherwise, set it to the maximum denominator that can
    /// represent the coefficients in rational form. This is a useful information to
    /// eliminate round-off errors when assembling the regularized coefficient matrix.
    unsigned max_denominator = 0;
};

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

    /// Set the options for regularizing linear constraints.
    auto setOptions(const RegularizerOptions& options) -> void;

    /// Regularize the optimum problem, state, and options before they are used in an optimization calculation.
    /// @param problem The optimum problem to be regularized.
    /// @param state The optimum state to be regularized.
    /// @param options The optimum options to be regularized.
    auto regularize(OptimumProblem& problem, OptimumState& state, OptimumOptions& options) -> void;

    /// Regularize the vectors `dg/dp` and `db/dp`, where `g = grad(f)`.
    auto regularize(Vector& dgdp, Vector& dbdp) -> void;

    /// Recover an optimum state to an state that corresponds to the original optimum problem.
    /// @param state[in,out] The optimum state regularized in method `regularize`.
    auto recover(OptimumState& state) -> void;

    /// Recover the sensitivity derivative `dxdp`.
    auto recover(Vector& dxdp) -> void;

  private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
