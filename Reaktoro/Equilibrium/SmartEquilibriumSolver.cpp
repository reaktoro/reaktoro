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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>

namespace Reaktoro {

struct SmartEquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The options for the smart equilibrium calculation
    SmartEquilibriumOptions options;

    /// The result of the last smart equilibrium calculation.
    SmartEquilibriumResult result;

    /// The solver for the equilibrium calculations
    EquilibriumSolver solver;

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

    /// A class used to store the node of tree for smart equilibrium calculations.
    struct TreeNode
    {
        Vector be;
        ChemicalState state;
        ChemicalProperties properties;
        EquilibriumSensitivity sensitivity;
    };

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::list<TreeNode> tree;

    /// The vector of amounts of species
    Vector n;

    /// Auxiliary vectors delta(n) and delta(lna)
    Vector dn, delta_lna;

    /// Construct a default SmartEquilibriumSolver::Impl instance.
    Impl()
    {}

    /// Construct an SmartEquilibriumSolver::Impl instance.
    Impl(const ChemicalSystem& system) : system(system), solver(system)
    {}

    /// Set the options for the equilibrium calculation.
    auto setOptions(const SmartEquilibriumOptions& options) -> void
    {
        this->options = options;
        solver.setOptions(options.learning);
    }

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition) -> void
    {
        solver.setPartition(partition);
    }

    // auto showTree(const Index& step) -> void
    // {
    //     /*
    //     std::cout << "\nTree: " << std::endl;
    //     for(auto node : tree_)
    //     {
    //         const auto be0_ = node.be;
    //         const ChemicalState state0_ = node.state;
    //         const ChemicalProperties properties0_ = node.properties;
    //         const EquilibriumSensitivity sensitivity0_ = node.sensitivity;
    //         const Index num_predicted = node.num_predicted;
    //         const Index step = node.step;
    //         const Index cell = node.cell;
    //         std::cout << "node ( step " << node.step
    //                   << ", cell " <<  node.cell << ") : "
    //                   << num_predicted << " predicted states" <<  std::endl;

    //     }
    //     std::cout <<  std::endl;
    //     */

    //     // Document times
    //     // ----------------------------------------------------------------------
    //     // ----------------------------------------------------------------------
    //     // The output stream of the data file
    //     std::ofstream datafile;
    //     std::string file = "../tree-" + std::to_string(step) + ".txt";
    //     //datafile.open("tree-" + std::to_string(step) + ".txt", std::ofstream::out | std::ofstream::trunc);
    //     datafile.open(file, std::ofstream::out | std::ofstream::trunc);
    //     // Output statuses such that steps' go vertically and cells horizontally
    //     if (datafile.is_open()) {

    //         // Output tree
    //         for(auto node : tree)
    //         {
    //             const auto be0_ = node.be;
    //             const ChemicalState state0_ = node.state;
    //             const ChemicalProperties properties0_ = node.properties;
    //             const EquilibriumSensitivity sensitivity0_ = node.sensitivity;
    //             const Index num_predicted = node.num_predicted;
    //             const Index step = node.step;
    //             const Index cell = node.cell;
    //             datafile << "node ( step " << node.step
    //                      << ", cell " <<  node.cell << ") : "
    //                      << num_predicted << " predicted states" <<  std::endl;

    //         }
    //         datafile <<  std::endl;
    //     }

    //     // Close the data file
    //     datafile.close();
    // }

    /// Learn how to perform a full equilibrium calculation (with tracking)
    auto learn(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        // Calculate the equilibrium state using conventional Gibbs energy minimization approach
        timeit( solver.solve(state, T, P, be),
            result.timing.learning_gibbs_energy_minimization= );

        // Store the result of the Gibbs energy minimization calculation performed during learning
        result.learning.gibbs_energy_minimization = solver.result();

        // Update the chemical properties of the system
        properties = solver.properties();

        // Store the computed solution into the knowledge tree
        timeit( tree.push_back({be, state, properties, solver.sensitivity()}),
            result.timing.learning_storage= );
    }

    /// Estimate the equilibrium state using sensitivity derivatives (profiling the expences)
    auto estimate(ChemicalState& state, double T, double P, VectorConstRef be) -> void
    {
        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

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

        // Relative and absolute tolerance parameters
        const auto reltol = options.reltol;
        const auto abstol = options.abstol;

        //---------------------------------------------------------------------------------------
        // Step 1: Search for the reference element (closest to the new state input conditions)
        //---------------------------------------------------------------------------------------
        tic(0);

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        {
            const auto& be_a = a.be;
            const auto& be_b = b.be;
            return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };

        // Find the entry with minimum "input" distance
        auto it = std::min_element(tree.begin(), tree.end(), distancefn);

        toc(0, result.timing.estimate_search);

        //----------------------------------------------------------------------------
        // Step 2: Calculate predicted state with a first-order Taylor approximation
        //----------------------------------------------------------------------------
        tic(1);

        // Get all the data stored in the reference element
        const auto& be0 = it->be;
        const auto& state0 = it->state;
        const auto& properties0 = it->properties;
        const auto& sensitivity0 = it->sensitivity;
        const auto& n0 = state0.speciesAmounts();

        // Get the sensitivity derivatives dln(a) / dn
        const auto lna_vec = properties0.lnActivities();
        const auto& dlnadn = lna_vec.ddn; // TODO: this line is assuming all species are equilibrium specie! get the rows and columns corresponding to equilibrium species
        const auto& lna0 = lna_vec.val;

        // Calculate perturbation of n
        dn.noalias() = sensitivity0.dndb * (be - be0); // delta(n) = dn/db * (b - b0)
        n.noalias() = n0 + dn;                         // n = n0 + delta(n)
        delta_lna.noalias() = dlnadn * dn;             // delta(ln(a)) = d(lna)/dn * delta(n)

        toc(1, result.timing.estimate_mat_vec_mul);

        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------
        tic(2);

        const auto& x = properties0.moleFractions();

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
        for(Index i = 0; i < n.size(); ++i)
        {
            // If the fraction is too small, skip the variational check
            if(x[i] < fraction_tol)
                continue;

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

        toc(2, result.timing.estimate_acceptance);


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

        // Check if smart estimation failed with respect to variation of chemical potentials or amounts
        if(!variation_check || !amount_check)
            return;

        // Set the output chemical state to the approximate amounts of species
        n.noalias() = abs(n); // TODO: this mirroring should be replaced by setting a very small mole amount for negative small amounts
        state.setSpeciesAmounts(n);

        // Set the estimate accepted status to true
        result.estimate.accepted = true;

        // Update the chemical properties of the system
        properties = properties0;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
    {
        tic(0);

        // Reset the result of the last smart equilibrium calculation
        result = {};

        // Perform a smart estimate of the chemical state
        timeit( estimate(state, T, P, be),
            result.timing.estimate= );

        // Perform a learning step if the smart prediction is not sactisfatory
        if(!result.estimate.accepted)
            timeit( learn(state, T, P, be), result.timing.learning= );

        toc(0, result.timing.solve);

        return result;
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

auto SmartEquilibriumSolver::setOptions(const SmartEquilibriumOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto SmartEquilibriumSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult
{
    return solve(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto SmartEquilibriumSolver::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto SmartEquilibriumSolver::result() const -> const SmartEquilibriumResult&
{
    return pimpl->result;
}

} // namespace Reaktoro

