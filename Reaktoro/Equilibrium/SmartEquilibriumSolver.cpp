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
#include <deque>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>

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

    /// The amounts of the elements in the equilibrium partition
    Vector be;

    /// The solution of the equilibrium problem
    Vector n, y, z, x, u;

    /// Auxiliary vectors
    Vector dn, dy, dz;

    /// A class used to store the node of tree for smart equilibrium calculations.
    struct TreeNode
    {
        Vector be;
        ChemicalState state;
        ChemicalProperties properties;
        EquilibriumSensitivity sensitivity;
    };

    /// The tree used to save the calculated equilibrium states and respective sensitivities
    std::deque<TreeNode> tree;

    /// Construct a default SmartEquilibriumSolver::Impl instance.
    Impl()
    {}

    /// Construct an SmartEquilibriumSolver::Impl instance with given partition.
    Impl(const Partition& partition_)
    : partition(partition_), system(partition_.system()), solver(partition_)
    {
    }

    /// Set the options for the equilibrium calculation.
    auto setOptions(const SmartEquilibriumOptions& options_) -> void
    {
        options = options_;

        options.learning.hessian = GibbsHessian::Exact; // Ensure the use of an exact Hessian of the Gibbs energy function

        solver.setOptions(options.learning);
    }

    /// Set the partition of the chemical system.
    auto setPartition(const Partition& partition_) -> void
    {
        partition = partition_;
        solver.setPartition(partition_);
    }

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
        result.estimate.accepted = false;

        // Skip estimation if no previous full computation has been done
        if(tree.empty())
            return;

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
        const auto& T0 = state0.temperature();
        const auto& P0 = state0.pressure();
        const auto& n0 = state0.speciesAmounts();
        const auto& y0 = state0.elementDualPotentials();
        const auto& z0 = state0.speciesDualPotentials();
        const auto u0 = properties0.chemicalPotentials();
        const auto x0 = properties0.moleFractions();



        const double RT = universalGasConstant*T0;


        // Calculate perturbation of n
        dn.noalias() = sensitivity0.dndb * (be - be0);
        dy.noalias() = sensitivity0.dydb * (be - be0);
        dz.noalias() = sensitivity0.dzdb * (be - be0);

        n.noalias() = n0 + dn;
        // y.noalias() = y0 + dy;
        // z.noalias() = z0 + dz;
        y.noalias() = y0;
        z.noalias() = z0;
        // y.noalias() = y0 + dy * RT;
        // z.noalias() = z0 + dz * RT;

        u = u0.val + u0.ddn * dn; // u = u_ref + dudn_ref * dn
        x = x0.val + x0.ddn * dn; // x = x_ref + dxdn_ref * dn

        // MatrixConstRef Ae = partition.formulaMatrixEquilibriumPartition();

        Matrix A = system.formulaMatrix();

        Vector r = abs(u - tr(A)*y - z);  // TODO: We should actually collect the entries in u and z corresponding to equilibrium species
        r /= universalGasConstant*T0;
        // Vector r = abs(sensitivity0.drdb * (be - be0));

        Index ispecies;

        const double error = (x.cwiseProduct(r)).maxCoeff(&ispecies);

        bool is_error_acceptable = error < 1e-1;

        // if(is_error_acceptable == false)
        {
            std::cout << "*****************************" << std::endl;
            std::cout << "*** Failed Smart Estimate ***" << std::endl;
            std::cout << "*****************************" << std::endl;
            std::cout << "Error = " << error << std::endl;
            std::cout << "Triggered by species = " << system.species(ispecies).name() << std::endl;
            std::cout << "-----------------------------" << std::endl;
            std::cout << std::left << std::setw(25) << "Species";
            std::cout << std::left << std::setw(25) << "r[i]";
            std::cout << std::left << std::setw(25) << "x[i]";
            std::cout << std::left << std::setw(25) << "r[i] * x[i]";
            std::cout << std::endl;
            std::cout << std::scientific;
            for(auto i = 0; i < r.size(); ++i)
            {
                if(i == ispecies) std::cout << "***************" << std::endl;
                std::cout << std::left << std::setw(25) << system.species(i).name();
                std::cout << std::left << std::setw(25) << r[i];
                std::cout << std::left << std::setw(25) << x[i];
                std::cout << std::left << std::setw(25) << r[i] * x[i];
                std::cout << std::endl;
                if(i == ispecies) std::cout << "***************" << std::endl;
            }
            std::cout << "=============================" << std::endl;
        }

        toc(1, result.timing.estimate_mat_vec_mul);

        //----------------------------------------------
        // Step 3: Checking the acceptance criterion
        //----------------------------------------------
        tic(2);

        // Correct negative mole numbers
        for(auto i = 0; i < n.size(); ++i)
            if(n[i] < 0.0) n[i] = 1e-12;  // TODO: 1e-12 should be replaced by a value provided via options.

        state.setSpeciesAmounts(n);
        state.setElementDualPotentials(y0 + dy * RT);
        state.setSpeciesDualPotentials(z0 + dz * RT);

        toc(2, result.timing.estimate_acceptance);

        // Check if smart estimation failed with respect to variation of chemical potentials or amounts
        // if(!is_error_acceptable || !amount_check)
        if(is_error_acceptable == false)
            return;

        // Set the estimate accepted status to true
        result.estimate.accepted = true;

        // Update the chemical properties of the system
        properties = properties0;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
    }

    /// Solve the equilibrium problem with given initial state
    auto solve(ChemicalState& state) -> SmartEquilibriumResult
    {
        const auto& iee = partition.indicesEquilibriumElements();
        const auto& ies = partition.indicesEquilibriumSpecies();
        const auto T = state.temperature();
        const auto P = state.pressure();
        be = state.elementAmountsInSpecies(ies)(iee);
        return solve(state, T, P, be);
    }

    /// Solve the equilibrium problem with given problem definition
    auto solve(ChemicalState& state, const EquilibriumProblem& problem) -> SmartEquilibriumResult
    {
        setPartition(problem.partition());
        const auto T = problem.temperature();
        const auto P = problem.pressure();
        const auto& iee = partition.indicesEquilibriumElements(); // This statement needs to be after setPartition
        be = problem.elementAmounts()(iee);
        return solve(state, T, P, be);
    }

    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> SmartEquilibriumResult
    {
        tic(0);

        // Absolutely ensure an exact Hessian of the Gibbs energy function is used in the calculations
        setOptions(options);

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
: pimpl(new Impl(Partition(system)))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(const Partition& partition)
: pimpl(new Impl(partition))
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
    return pimpl->solve(state, problem);
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

