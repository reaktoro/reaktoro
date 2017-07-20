// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "SmartEquilibriumSolver.hpp"

// C++ includes
#include <list>
#include <tuple>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktoro {

struct SmartEquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::list<std::tuple<Vector, ChemicalState, EquilibriumSensitivity>> tree;

    /// Construct a default SmartEquilibriumSolver::Impl instance.
    Impl()
    {}

    /// Construct an SmartEquilibriumSolver::Impl instance.
    Impl(const ChemicalSystem& system)
    : system(system), solver(system)
    {}

    /// Set the options for the equilibrium calculation.
    auto setOptions(const EquilibriumOptions& options) -> void
    {
        this->options = options;
        solver.setOptions(options);
    }

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition) -> void
    {
        solver.setPartition(partition);
    }

    /// Learn how to perform a full equilibrium calculation.
    auto learn(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
    {
        EquilibriumResult res = solver.solve(state, T, P, be);
        tree.emplace_back(be, state, solver.sensitivity());
        return res;
    }

    auto estimate(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
    {
        if(tree.empty())
            return {};

        using TreeNodeType = std::tuple<Vector, ChemicalState, EquilibriumSensitivity>;

        EquilibriumResult res;

        auto comp = [&](const TreeNodeType& a, const TreeNodeType& b)
        {
            const Vector& be_a = std::get<0>(a);
            const Vector& be_b = std::get<0>(b);
            return norm(be_a - be) < norm(be_b - be);
        };

        auto it = std::min_element(tree.begin(), tree.end(), comp);

        const Vector& be0 = std::get<0>(*it);
        const ChemicalState& state0 = std::get<1>(*it);
        const EquilibriumSensitivity& sensitivity0 = std::get<2>(*it);
        const Vector& n0 = state0.speciesAmounts();
        Vector n = n0 + sensitivity0.dnedbe * (be - be0);

        const auto reltol = options.smart.reltol;
        const auto abstol = options.smart.abstol;

        if(((n - n0).array().abs() <= abstol + reltol*n0.array().abs()).all())
        {
            state.setSpeciesAmounts(n);
            res.optimum.succeeded = true;
            res.smart = true;
            return res;
        }

        return res;
    }

    auto solve(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
    {
        EquilibriumResult res = estimate(state, T, P, be);

        if(res.optimum.succeeded) return res;

        res += learn(state, T, P, be);

        return res;
    }
};

SmartEquilibriumSolver::SmartEquilibriumSolver()
: pimpl(new Impl())
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(const SmartEquilibriumSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

auto SmartEquilibriumSolver::operator=(SmartEquilibriumSolver other) -> SmartEquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

SmartEquilibriumSolver::~SmartEquilibriumSolver()
{}

auto SmartEquilibriumSolver::setOptions(const EquilibriumOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto SmartEquilibriumSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto SmartEquilibriumSolver::learn(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
{
    return pimpl->learn(state, T, P, be);
}

auto SmartEquilibriumSolver::learn(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return learn(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto SmartEquilibriumSolver::estimate(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
{
    return pimpl->estimate(state, T, P, be);
}

auto SmartEquilibriumSolver::estimate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return estimate(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, const Vector& be) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return solve(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

} // namespace Reaktoro

