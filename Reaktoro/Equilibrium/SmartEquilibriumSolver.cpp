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
#include <algorithm>
#include <deque>
#include <tuple>

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/LU>

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
#include <Reaktoro/Optimization/Canonicalizer.hpp>

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
    Vector n, y, z, x, u, r;

    /// Auxiliary vectors
    Vector dn, dy, dz;

    std::deque<Index> priority;

    std::deque<Index> ranking;





    Vector a_dn;
    Vector a_dy;
    Vector a_dz;
    Vector a_n0;
    Vector a_y0;
    Vector a_z0;
    Vector a_n;
    Vector a_y;
    Vector a_z;
    Vector b_dn;
    Vector b_dy;
    Vector b_dz;
    Vector b_n0;
    Vector b_y0;
    Vector b_z0;
    Vector b_n;
    Vector b_y;
    Vector b_z;
    Vector a_u;
    Vector b_u;
    Vector a_x;
    Vector b_x;
    Vector a_r;
    Vector b_r;

    /// A class used to store the node of tree for smart equilibrium calculations.
    struct TreeNode
    {
        Vector be;
        Vector bebar;
        ChemicalState state;
        ChemicalProperties properties;
        EquilibriumSensitivity sensitivity;
        Matrix Mb;
        VectorXi imajor;
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

        // Tweak the options for the Gibbs energy minimization during learning operations.
        options.learning.hessian = GibbsHessian::Exact; // ensure the use of an exact Hessian of the Gibbs energy function
        options.learning.optimum.tolerance = 1e-10; // ensure the use of a stricter residual tolerance for the Gibbs energy minimization

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
        priority.push_back(priority.size());

        ranking.push_back(0);

        const auto RT = universalGasConstant * T;

        const auto& n = state.speciesAmounts();

        const auto u = properties.chemicalPotentials();
        const auto x = properties.moleFractions();

        const auto& dndb = solver.sensitivity().dndb;
        const auto& dudn = u.ddn;

        const Matrix dudb = dudn * dndb;

        const double nsum = sum(n);

        const auto& A = system.formulaMatrix();

        Canonicalizer canonicalizer(A);
        canonicalizer.updateWithPriorityWeights(n);

        const auto eps_n = options.amount_fraction_cutoff * nsum;
        const auto eps_x = options.mole_fraction_cutoff;

        const auto S = canonicalizer.S();

        auto imajorminor = canonicalizer.Q();
        const auto nummajor = std::partition(imajorminor.begin(), imajorminor.end(),
            [&](Index i) { return n[i] >= eps_n && x.val[i] >= eps_x; }) - imajorminor.begin();

        const auto imajor = imajorminor.head(nummajor);

        // const auto iprimary = canonicalizer.indicesBasicVariables();
        // const auto isecondary = canonicalizer.indicesNonBasicVariables();

        // Indices imajor;
        // imajor.reserve(isecondary.size());
        // for(auto i = 0; i < isecondary.size(); ++i)
        //     if(n[isecondary[i]] >= eps_n && x.val[isecondary[i]] >= eps_x)
        //         imajor.push_back(i);

        // VectorXi isecondary_major = isecondary(imajor);

        const Vector um = u.val(imajor);
        const auto dumdb = rows(dudb, imajor);

        const Matrix Mb = diag(inv(um)) * dumdb;

        timeit( tree.push_back({be, be/sum(be), state, properties, solver.sensitivity(), Mb, imajor}),
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

        Vector rs;
        Vector n, ns;
        Vector du, dus, dns;
        Vector x;
        Vector bebar = be/sum(be);
        Vector dbe;

        double nmin, nsum, uerror;
        Index inmin, iuerror;

        double nerror;
        Index inerror;


        // Check if an entry in the database pass the error test.
        // It returns (`success`, `error`, `ispecies`), where
        //   - `success` is true if error test succeeds, false otherwise.
        //   - `error` is the first error violating the tolerance
        //   - `ispecies` is the index of the species that fails the error test
        auto pass_error_test = [&](const auto& node) -> std::tuple<bool, double, Index>
        {
            using std::abs;
            using std::max;
            const auto& be0 = node.be;
            const auto& Mb0 = node.Mb;
            const auto& imajor = node.imajor;

            dbe.noalias() = be - be0;

            double error = 0.0;
            for(auto i = 0; i < imajor.size(); ++i) {
                error = max(error, abs(Mb0.row(i) * dbe));
                if(error >= reltol)
                    return { false, error, imajor[i] };
            }

            return { true, error, n.size() };
        };

        for(auto inode : priority)
        {
            const auto& node = tree[inode];
            const auto& imajor = node.imajor;

            const auto [success, error, ispecies] = pass_error_test(node);

            if(success)
            {
                const auto& be0 = node.be;
                const auto& n0 = node.state.speciesAmounts();
                const auto& y0 = node.state.elementDualPotentials();
                const auto& z0 = node.state.speciesDualPotentials();
                const auto& dndb0 = node.sensitivity.dndb;

                n.noalias() = n0 + dndb0 * (be - be0);

                nmin = min(n(imajor));
                nsum = sum(n);

                const auto eps_n = options.amount_fraction_cutoff * nsum;

                if(nmin < -eps_n)
                    continue;

                ranking[inode] += 1;

                std::sort(priority.begin(), priority.begin() + inode + 1,
                    [&](Index l, Index r) { return ranking[l] > ranking[r]; });

                state.setSpeciesAmounts(n);
                state.setElementDualPotentials(y0);
                state.setSpeciesDualPotentials(z0);

                // Update the chemical properties of the system
                properties = node.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
                result.estimate.accepted = true;
                return;
            }
            else continue;
        }

        result.estimate.accepted = false;
        return;
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
            timeit( learn(state, T, P, be), result.timing.learn= );

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

