// Reaktoro is a unified framework for modeling chemically reactive systems.
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
#include <Reaktoro/Common/Profiling.hpp>
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
    std::list<TreeNode> tree;

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
            for(auto node : tree)
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

        // Profiling time variables
        profiling( Time start_learn; );
        profiling( Time start_parts; );

        // Start profiling learning
        profiling( start_learn = time(); );
        profiling( start_parts = time(); );

        // Construct chemical state by the conventional approach
        res = solver.solve(state, T, P, be);

        // Profiling gibbs minimization
        profiling( res.smart.learn_stats.time_gibbs_min = elapsed(start_parts); );

        profiling( start_parts = time(); );

        // Add the reference state into the storage
        tree.emplace_back(be, state, solver.properties(), solver.sensitivity(), 0, step, cell);

        // Profiling ref. element store
        profiling( res.smart.learn_stats.time_store = elapsed(start_parts); );

        // Stop profiling learning
        profiling( res.smart.learn_stats.time_learn = elapsed(start_learn); );

        return res;
    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
    {
        // If the tree is empty skip estimation
        if(tree.empty())
            return {};

        // Class that stores info about the equilibrium computations
        EquilibriumResult res;

        // Profiling time variables
        profiling( Time start_estimate; );
        profiling( Time start_parts; );

        // Start profiling estimating
        profiling( start_estimate = time(); );

        // Comparison function based on the Euclidean distance
        auto comp_ = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& be_a = a.be;
            const auto& be_b = b.be;
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();
        };

        // Step 1: search for the reference element (closest to the new state be)

        // Profiling ref. element search
        profiling( start_parts = time(); );

        // Find the reference element (closest to the new state be)
        auto it = std::min_element(tree.begin(), tree.end(), comp_);

        profiling( res.smart.estimate_stats.time_search = elapsed(start_parts); );

        // Get all the data stored in the reference element
        const auto& be0 = it->be;
        const ChemicalState& state0 = it->state;
        const ChemicalProperties& properties0 = it->properties;
        const EquilibriumSensitivity& sensitivity0 = it->sensitivity;

        // const auto& n0 = state0.speciesAmounts();
        const auto& n0 = state0.speciesAmounts();

        // Get the sensitivity derivatives dln(a) / dn
        MatrixConstRef dlnadn = properties0.lnActivities().ddn; // TODO: this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
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

        // Step 2: calculate predicted state

        // Profiling matrix-vector manipulations
        profiling( start_parts = time(); );

        // Calculate perturbation of n
        dn.noalias() = sensitivity0.dndb * (be - be0);    // delta(n) = dn/db * (b - b0)
        n.noalias() = n0 + dn;                              // n = n0 + delta(n)
        delta_lna.noalias() = dlnadn * dn;                  // delta(ln(a)) = d(lna)/dn * delta(n)

        const auto& x = properties0.moleFractions();

        profiling( res.smart.estimate_stats.time_mat_vect_mult = elapsed(start_parts); );

        // Step 3: checking the acceptance criterion

        // Profiling acceptance test
        profiling( start_parts = time(); );

        // The estimated activity of all species must be positive.
        // This results in the need to check delta(ln(a[i])) > -1 for all species.

        // Perform the check for the negative amounts
        const double cutoff = -1e-5;    // relaxation parameter
        const bool amount_check = n.minCoeff() > cutoff;

        // Assign small values to all the amount in the interval [cutoff, 0]
        // (instead of mirroring bellow)
        // for(int i = 0; i < n.size(); ++i) if(n[i] > cutoff && n[i] < 0) n[i] = options.epsilon;

        // The estimated ln(a[i]) of each species must not be
        // too far away from the reference value ln(aref[i])
        // const bool variation_check = (delta_lna.array().abs() <= abstol + reltol * lna0.array().abs()).all();
        ///*
        // Check in the loop mole fractions and variations of the ln(a)
        ///*
        bool variation_check = true;
        const double fraction_tol = abstol * 1e-2;
        for(int i = 0; i < n.size(); ++i)
        {
            // If the fraction is too small, skip the variational check
            if(x[i] < fraction_tol) continue;

            // Perform the variational check
            if(std::abs(delta_lna[i]) > abstol + reltol * std::abs(lna0[i])) {
                variation_check = false;    // variation test failed
                /*
                // Info about the failed species
                std::cout << "i: " << i << std::endl;
                std::cout << "Species: " << system.species(i).name() << std::endl;
                std::cout << "Amount(new): " << n[i] << std::endl;
                std::cout << "Amount(old): " << n0[i] << std::endl;
                std::cout << "logActivity(new): " << lna0[i] + delta_lna[i] << std::endl;
                std::cout << "logActivity(old): " << lna0[i] << std::endl;
                std::cout << "RelativeDifference(Amount): " << std::abs((n[i] - n0[i])/n0[i]) << std::endl;
                std::cout << "RelativeDifference(lnActivity): " << std::abs(delta_lna[i])/(std::abs(lna0[i]) + 1.0) << std::endl;
                std::cout << "VariationCheck: " << variation_check << std::endl;
                std::cout << "AmountCheck: " << amount_check << std::endl;
                */
            }
        }
        //*/
        //*/
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
            n.noalias() = abs(n); // this mirroring is replaced by setting a very small mole amount
            state.setSpeciesAmounts(n);
            res.optimum.succeeded = true;
            res.smart.succeeded = true;
            it->num_predicted ++;
        }

        profiling( res.smart.estimate_stats.time_acceptance = elapsed(start_parts); );

        // Stop profiling learning
        profiling( res.smart.estimate_stats.time_estimate = elapsed(start_estimate); );

        return res;
    }

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
auto SmartEquilibriumSolver::showTree(const Index & step) const -> void
{
    return pimpl->showTree(step);
}

} // namespace Reaktoro

