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

#include "SmartEquilibriumSolver.hpp"

// C++ includes
#include <iostream> // todo remove
#include <list>
#include <tuple>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
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
    std::list<std::tuple<Vector, ChemicalState, ChemicalProperties, EquilibriumSensitivity>> tree;

    /// The vector of amounts of species
    Vector n;

    Vector dn, delta_lna;

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
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
    {
        EquilibriumResult res = solver.solve(state, T, P, be);
        tree.emplace_back(be, state, solver.properties(), solver.sensitivity());
        return res;
    }

    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
    {
        if(tree.empty())
            return {};

        using TreeNodeType = std::tuple<Vector, ChemicalState, ChemicalProperties, EquilibriumSensitivity>;

        EquilibriumResult res;

        auto comp = [&](const TreeNodeType& a, const TreeNodeType& b) {
            const auto& be_a = std::get<0>(a);
            const auto& be_b = std::get<0>(b);
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();
        };

        auto it = std::min_element(tree.begin(), tree.end(), comp);

        const auto& be0 = std::get<0>(*it);
        const ChemicalState& state0 = std::get<1>(*it);
        const ChemicalProperties& properties0 = std::get<2>(*it);
        const EquilibriumSensitivity& sensitivity0 = std::get<3>(*it);
        const auto& n0 = state0.speciesAmounts();

        MatrixConstRef dlnadn = properties0.lnActivities().ddn; // TODO this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
        const auto& lna0 = properties0.lnActivities().val;

        // TODO Fixing negative amounts
        // Once some species are found to have negative values, first check
        // After the projection, assume species i has negative amounts.
        // 1) Check if n(i,new) is greater than, say, -1.0e-4.
        //    If so, just approximate n(i,new) to, say, 1e-25 (some small number
        //    below abstol!).
        // 2)

        //        const double ln10 = 2.302585;

        //        const auto reltol = ln10 * options.smart.reltol;
        //        const auto abstol = ln10 * options.smart.abstol;
        const auto reltol = options.smart.reltol;
        const auto abstol = options.smart.abstol;

        //        n = n0 + sensitivity0.dnedbe * (be - be0);
        dn.noalias() = sensitivity0.dndb * (be - be0); // n is actually delta(n)

        n.noalias() = n0 + dn;

        delta_lna.noalias() = dlnadn * dn;

        // The estimated ln(a[i]) of each species must not be
        // too far away from the reference value ln(aref[i])
        const bool variation_check = (delta_lna.array().abs() <=
                                      abstol + reltol * lna0.array().abs())
                                         .all();

        // The estimated activity of all species must be positive.
        // This results in the need to check delta(ln(a[i])) > -1 for all species.
        //        const bool activity_check = delta_lna.minCoeff() > -1.0;
        const bool amount_check = n.minCoeff() > -1e-5;
        //
        //        for(int i = 0; i < n.size(); ++i)
        //        {
        //            if(n[i] > 0.0) continue;
        //            if(n[i] > -1e-5) n[i] = abstol;
        //        }

        //        if(variation_check && activity_check && amount_check)
        if(variation_check && amount_check)
        //        if(variation_check)
        //        if(((n - n0).array().abs() <= abstol + reltol*n0.array().abs()).all())
        {
            n.noalias() = abs(n); // TODO abs needs only to be applied to negative values
            state.setSpeciesAmounts(n);
            res.optimum.succeeded = true;
            res.smart.succeeded = true;
            return res;
        }

        // std::cout << "=======================" << std::endl;
        // std::cout << "Smart Estimation Failed" << std::endl;
        // std::cout << "=======================" << std::endl;
        // for(int i = 0; i < n.size(); ++i)
        // {
        //     const bool variation_check = (std::abs(delta_lna[i]) <=
        //             abstol + reltol * std::abs(lna0[i]));
        //     const bool activity_check = delta_lna[i] > -1.0;
        //     const bool amount_check = n[i] > -1e-5;

        //     if(variation_check && activity_check && amount_check)
        //         continue;

        //     std::cout << "Species: " << system.species(i).name() << std::endl;
        //     std::cout << "Amount(new): " << n[i] << std::endl;
        //     std::cout << "Amount(old): " << n0[i] << std::endl;
        //     std::cout << "logActivity(new): " << lna0[i] + delta_lna[i] << std::endl;
        //     std::cout << "logActivity(old): " << lna0[i] << std::endl;
        //     std::cout << "RelativeDifference(Amount): " << std::abs((n[i] - n0[i])/n0[i]) << std::endl;
        //     std::cout << "RelativeDifference(lnActivity): " << std::abs(delta_lna[i])/(std::abs(lna0[i]) + 1.0) << std::endl;
        //     std::cout << "VariationCheck: " << variation_check << std::endl;
        //     std::cout << "ActivityCheck: " << activity_check << std::endl;
        //     std::cout << "AmountCheck: " << amount_check << std::endl;
        //     std::cout << std::endl;
        // }

        return res;
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
    {
        EquilibriumResult res = estimate(state, T, P, be);

        if(res.optimum.succeeded)
            return res;

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

auto SmartEquilibriumSolver::learn(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
{
    return pimpl->learn(state, T, P, be);
}

auto SmartEquilibriumSolver::learn(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return learn(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto SmartEquilibriumSolver::estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
{
    return pimpl->estimate(state, T, P, be);
}

auto SmartEquilibriumSolver::estimate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return estimate(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return solve(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto SmartEquilibriumSolver::properties() const -> const ChemicalProperties&
{
    RuntimeError("Could not calculate the chemical properties.",
                 "This method has not been implemented yet.");
}

} // namespace Reaktoro
