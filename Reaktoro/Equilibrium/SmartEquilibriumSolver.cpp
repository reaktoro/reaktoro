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
#include <list>
#include <fstream>
#include <numeric>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/TimeUtils.hpp>

#include <Reaktoro/Core/ChemicalProperties.hpp>
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
    //std::list<std::tuple<Vector, ChemicalState, ChemicalProperties, EquilibriumSensitivity>> tree;

    /// A class used to store the node of tree for smart equilibrium calculations.
    struct TreeNode{

        Vector be;
        ChemicalState state;
        ChemicalProperties properties;
        EquilibriumSensitivity sensitivity;

        Index num_predicted;
        Index step;
        Index cell;

        TreeNode(Vector be,
                 ChemicalState state,
                 ChemicalProperties properties,
                 EquilibriumSensitivity sensitivity,
                 Index num_predicted, Index step, Index cell)
                : be(be), state(state), properties(properties), sensitivity(sensitivity),
                  num_predicted(num_predicted), step(step), cell(cell){}
    };

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::list<TreeNode> tree_;

    /// The vector of amounts of species
    Vector n;

    Vector dn, delta_lna;

    /// Construct a default SmartEquilibriumSolver::Impl instance.
    Impl()
    {}

    /// Construct an SmartEquilibriumSolver::Impl instance.
    Impl(const ChemicalSystem& system) : system(system), solver(system)
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

    auto showTree(const Index& step) -> void
    {
        /*
        std::cout << "\nTree: " << std::endl;
        for(auto node : tree_)
        {
            const auto be0_ = node.be;
            const ChemicalState state0_ = node.state;
            const ChemicalProperties properties0_ = node.properties;
            const EquilibriumSensitivity sensitivity0_ = node.sensitivity;
            const Index num_predicted = node.num_predicted;
            const Index step = node.step;
            const Index cell = node.cell;
            std::cout << "node ( step " << node.step
                      << ", cell " <<  node.cell << ") : "
                      << num_predicted << " predicted states" <<  std::endl;

        }
        std::cout <<  std::endl;
        */

        // Document times
        // ----------------------------------------------------------------------
        // ----------------------------------------------------------------------
        // The output stream of the data file
        std::ofstream datafile;
        std::string file = "../tree-" + std::to_string(step) + ".txt";
        //datafile.open("tree-" + std::to_string(step) + ".txt", std::ofstream::out | std::ofstream::trunc);
        datafile.open(file, std::ofstream::out | std::ofstream::trunc);
        // Output statuses such that steps' go vertically and cells horizontally
        if (datafile.is_open()) {

            // Output tree
            for(auto node : tree_)
            {
                const auto be0_ = node.be;
                const ChemicalState state0_ = node.state;
                const ChemicalProperties properties0_ = node.properties;
                const EquilibriumSensitivity sensitivity0_ = node.sensitivity;
                const Index num_predicted = node.num_predicted;
                const Index step = node.step;
                const Index cell = node.cell;
                datafile << "node ( step " << node.step
                         << ", cell " <<  node.cell << ") : "
                         << num_predicted << " predicted states" <<  std::endl;

            }
            datafile <<  std::endl;
        }

        // Close the data file
        datafile.close();
    }

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be, Index step, Index cell) -> EquilibriumResult
    {
        EquilibriumResult res;
        Time start_learn, start_parts;

        // Start profiling learning
        if (options.track_statistics) start_learn = time();

        if (options.track_statistics) start_parts = time();
        // Construct chemical state by the conventional approach
        res = solver.solve(state, T, P, be);
        if (options.track_statistics)
            res.smart.learn_stats.time_gibbs_min = elapsed(start_parts); // Profiling gibbs minimization

        if (options.track_statistics) start_parts = time();
        // Add the reference state into the storage
        //tree.emplace_back(be, state, solver.properties(), solver.sensitivity());
        //TreeNode new_state(be, state, solver.properties(), solver.sensitivity(), 0, step, cell);
        tree_.emplace_back(be, state, solver.properties(), solver.sensitivity(), 0, step, cell);
        if (options.track_statistics)
            res.smart.learn_stats.time_store = elapsed(start_parts); // Profiling ref. element store

        // Stop profiling learning
        if (options.track_statistics) res.smart.learn_stats.time_learn = elapsed(start_learn);

        return res;
    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
    {

        // If the tree is empty abort estimation
        //if(tree.empty())
        //    return {};
        if(tree_.empty())
            return {};

        //using TreeNodeType = std::tuple<Vector, ChemicalState, ChemicalProperties, EquilibriumSensitivity>;


        // Class that stores info about the equilibrium computations
        EquilibriumResult res;
        Time start_estimate, start_parts;

        // Start profiling estimating
        if (this->options.track_statistics) start_estimate = time();


        // Comparison function based on the Euclidean distance
        /*
        auto comp = [&](const TreeNodeType& a, const TreeNodeType& b)
        {
            const auto& be_a = std::get<0>(a);
            const auto& be_b = std::get<0>(b);
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();
        };
        */
        auto comp_ = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& be_a = a.be;
            const auto& be_b = b.be;
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();
        };

        // Step 1: search for the reference element (closest to the new state be)

        // Profiling ref. element search
        if (this->options.track_statistics) start_parts = time();

        // Find the reference element (closest to the new state be)
        //auto it = std::min_element(tree.begin(), tree.end(), comp);
        auto it_ = std::min_element(tree_.begin(), tree_.end(), comp_);

        if (this->options.track_statistics)
            res.smart.estimate_stats.time_search = elapsed(start_parts);

        // Get all the data stored in the reference element
        /*
        const auto& be0 = std::get<0>(*it);
        const ChemicalState& state0 = std::get<1>(*it);
        const ChemicalProperties& properties0 = std::get<2>(*it);
        const EquilibriumSensitivity& sensitivity0 = std::get<3>(*it);
        */
        const auto& be0_ = it_->be;
        const ChemicalState& state0_ = it_->state;
        const ChemicalProperties& properties0_ = it_->properties;
        const EquilibriumSensitivity& sensitivity0_ = it_->sensitivity;
        const Index& num_predicted = it_->num_predicted;

        // const auto& n0 = state0.speciesAmounts();
        const auto& n0 = state0_.speciesAmounts();

        // Get the sensitivity derivatives dln(a) / dn
        // MatrixConstRef dlnadn = properties0.lnActivities().ddn; // TODO: this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
        // const auto& lna0 = properties0.lnActivities().val;

        MatrixConstRef dlnadn = properties0_.lnActivities().ddn; // TODO: this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
        const auto& lna0 = properties0_.lnActivities().val;

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

        // Step 2: calculate predicted state

        // Profiling matrix-vector manipulations
        if (this->options.track_statistics) start_parts = time();

        // Calculate perturbation of n
        // n = n0 + sensitivity0.dnedbe * (be - be0);
        //dn.noalias() = sensitivity0.dndb * (be - be0); // n is actually delta(n)
        dn.noalias() = sensitivity0_.dndb * (be - be0_); // n is actually delta(n)
        n.noalias() = n0 + dn;
        delta_lna.noalias() = dlnadn * dn;

        const auto total_amount = sum(n);
        const auto mole_fractions = n / total_amount;

        /*
        std::cout << "total: " << total_amount << std::endl;
        std::cout << "mole fraction: \n " << mole_fractions << std::endl;
        std::cout << "amounts: \n " << n << std::endl;
        */

        if (this->options.track_statistics)
            res.smart.estimate_stats.time_mat_vect_mult = elapsed(start_parts);

        // Step 3: checking the acceptance criterion

        // Profiling acceptance test
        if (this->options.track_statistics) start_parts = time();

        // The estimated ln(a[i]) of each species must not be
        // too far away from the reference value ln(aref[i])
        const bool variation_check = (delta_lna.array().abs() <= abstol + reltol * lna0.array().abs()).all();
        //const auto variation_check_vector = delta_lna.array().abs() <= abstol + reltol * lna0.array().abs();
        //std::cout << variation_check_vector << std::endl;

        // The estimated activity of all species must be positive.
        // This results in the need to check delta(ln(a[i])) > -1 for all species.
        const bool activity_check = delta_lna.minCoeff() > -1.0;
        const bool amount_check = n.minCoeff() > -1e-5;
        const bool moll_fraction_check = mole_fractions.minCoeff() > -1e-5;
        //
        //        for(int i = 0; i < n.size(); ++i)
        //        {
        //            if(n[i] > 0.0) continue;
        //            if(n[i] > -1e-5) n[i] = abstol;
        //        }

        /*
        std::cout << "=======================" << std::endl;
        std::cout << "Smart Estimation Failed" << std::endl;
        std::cout << "=======================" << std::endl;
        for(int i = 0; i < n.size(); ++i)
        {
            const bool variation_check = (std::abs(delta_lna[i]) <= abstol + reltol * std::abs(lna0[i]));
            const bool activity_check = delta_lna[i] > -1.0;
            const bool amount_check = n[i] > -1e-5;
            const bool moll_fraction_check = n[i] / total_amount > -1e-5;

            if(variation_check && activity_check && amount_check)
                continue;

            std::cout << "i: " << i << std::endl;
            std::cout << "Species: " << system.species(i).name() << std::endl;
            std::cout << "Amount(new): " << n[i] << std::endl;
            std::cout << "Amount(old): " << n0[i] << std::endl;
            std::cout << "logActivity(new): " << lna0[i] + delta_lna[i] << std::endl;
            std::cout << "logActivity(old): " << lna0[i] << std::endl;
            std::cout << "RelativeDifference(Amount): " << std::abs((n[i] - n0[i])/n0[i]) << std::endl;
            std::cout << "RelativeDifference(lnActivity): " << std::abs(delta_lna[i])/(std::abs(lna0[i]) + 1.0) << std::endl;
            std::cout << "VariationCheck: " << variation_check << std::endl;
            std::cout << "ActivityCheck: " << activity_check << std::endl;
            std::cout << "AmountCheck: " << amount_check << std::endl;
            std::cout << "Mole fraction: " << moll_fraction_check << std::endl;

            std::cout << std::endl;
        }

        */

        if(variation_check && amount_check)
        {
            n.noalias() = abs(n);
            state.setSpeciesAmounts(n);
            res.optimum.succeeded = true;
            res.smart.succeeded = true;
            it_->num_predicted ++;
        }

        if (this->options.track_statistics)
            res.smart.estimate_stats.time_acceptance = elapsed(start_parts);

        // Stop profiling learning
        if (this->options.track_statistics)
            res.smart.estimate_stats.time_estimate = elapsed(start_estimate);

        return res;
    }
    /*
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
    {
        // Attempt to estimate the result by on-demand learning
        EquilibriumResult res = estimate(state, T, P, be);

        // If the obtained result satisfies the accuracy criterion, we accept it
        if(res.optimum.succeeded) return res;

        // Otherwise, trigger learning (conventional approach)
        res += learn(state, T, P, be);

        showTree();
        return res;
    }
    */
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be, Index step, Index cell) -> EquilibriumResult
    {
        // Attempt to estimate the result by on-demand learning
        EquilibriumResult res = estimate(state, T, P, be);

        // If the obtained result satisfies the accuracy criterion, we accept it
        if(res.optimum.succeeded) return res;

        // Otherwise, trigger learning (conventional approach)
        res += learn(state, T, P, be, step, cell);

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

auto SmartEquilibriumSolver::learn(ChemicalState& state, double T, double P, VectorConstRef be, Index step, Index cell) -> EquilibriumResult
{
    return pimpl->learn(state, T, P, be, step, cell);
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

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be, Index step, Index cell) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be, step, cell);
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
auto SmartEquilibriumSolver::showTree(const Index & step) const -> void{
    return pimpl->showTree(step);
}

} // namespace Reaktoro

