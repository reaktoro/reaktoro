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
        Indices isensitive;
        Indices imajor;
        ChemicalVector u;
        ChemicalVector x;
        Matrix dudb;
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

        const auto& R = universalGasConstant;

        const auto& n = state.speciesAmounts();

        const auto u = properties.chemicalPotentials();
        const auto x = properties.moleFractions();

        const auto& dndb = solver.sensitivity().dndb;
        const auto& dudn = u.ddn;

        const Matrix dudb = dudn * dndb;

        Matrix dnudb = diag(n) * dudb;
        dnudb += diag(u.val) * dndb;
        dnudb /= R*T;

        const double nsum = sum(n);

        const auto& A = system.formulaMatrix();


        // const Matrix W = diag(n/sum(n)) * abs(dudb);
        // Vector w = log(abs(n));
        // const auto wmin = min(w);
        // const auto wmax = max(w);
        // w.noalias() = (wmax - w)/(w - wmin);
        // w.noalias() = (w + 100.0)/(w + 1.0);
        // for(unsigned i = 0; i < w.size(); ++i)
        //     if(!std::isfinite(w[i])) w[i] = 1.0;

        // // Vector w = log(x.val);
        // // Vector w = log(n/nsum);
        // Vector w = log10(n/nsum);
        // const auto wmin = min(w);
        // const auto wmax = max(w);
        // w.noalias() = w + (wmax - wmin) + 1.0;

        // w.noalias() = (wmax - w)/(w - wmin);
        // w.noalias() = (w + 100.0)/(w + 1.0);
        // for(unsigned i = 0; i < w.size(); ++i)
        //     if(!std::isfinite(w[i])) w[i] = 1.0;


        // const Matrix W = diag(w) * A.transpose();

        // VectorXi ordering = VectorXi::LinSpaced(n.size(), 0, n.size());
        // std::sort(ordering.begin(), ordering.end(), [&](Index l, Index r) { return n[l] > n[r]; });

        // Matrix W(A.cols(), A.rows());
        // for(auto i = 0; i < W.rows(); ++i)
        //     W.row(i) = A.col(ordering[i]);

        // const Matrix W = diag(w) * dudb;
        // const Matrix W = diag(x.val) * abs(dudb);
        // const Matrix W = dudb;

        // const auto lu = W.fullPivLu();

        // const auto isensitive = lu.permutationP().indices().head(lu.rank());


        Canonicalizer canonicalizer(A);
        canonicalizer.updateWithPriorityWeights(n);

        const double epsmajor = 1e-14 * nsum;
        Indices imajor;
        imajor.reserve(n.size());
        for(auto i : canonicalizer.Q())
            if(n[i] >= epsmajor)
                imajor.push_back(i);

        Indices isensitive;
        isensitive.reserve(canonicalizer.numBasicVariables());
        for(auto i : canonicalizer.indicesBasicVariables())
            if(x.val[i] > options.mole_fraction_cutoff)
                isensitive.push_back(i);

        timeit( tree.push_back({be, be/sum(be), state, properties, solver.sensitivity(), isensitive, imajor, u, x, dudb}),
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

        Vector n, ns;
        Vector du, dus, dns;
        Vector x;
        Vector bebar = be/sum(be);

        double nmin, ntot, uerror;
        Index inmin, iuerror;

        double nerror;
        Index inerror;

        for(auto inode : priority)
        {
            const auto& node = tree[inode];
            const auto& n0 = node.state.speciesAmounts();
            const auto& be0 = node.be;
            const auto& bebar0 = node.bebar;
            const auto& u0 = node.u;
            const auto& dndb = node.sensitivity.dndb;
            const auto& imajor = node.imajor;
            const auto& is = node.isensitive;
            const auto& dudb0 = node.dudb;

            const auto us0 = u0.val(is);
            const auto dusdb0 = rows(dudb0, is);

            dus = dusdb0 * (be - be0);
            dus.noalias() = abs(dus / us0);

            // du = dudb0 * (be - be0);
            // du.noalias() = abs(du / u0.val);

            // std::cout << "------------------------------------------------------" << std::endl;
            // std::cout << "Sensitive Species:" << std::endl;
            // auto j = 0;
            // for(auto i : is) {
            //     auto name = system.species(i).name();
            //     std::cout << std::left << std::setw(15) << name;
            //     std::cout << std::left << std::setw(15) << dus[j];
            //     std::cout << std::left << std::setw(15) << du[i];
            //     std::cout << std::left << std::setw(15) << n0[i];
            //     std::cout << std::left << std::setw(15) << std::abs(u0.val[i])/(universalGasConstant * T);
            //     std::cout << std::endl;
            //     ++j;
            // }
            // std::cout << "------------------------------------------------------" << std::endl;

            // const auto ns0 = n0(is);
            // const auto us0 = u0.val(is);
            // const auto dnsdb0 = rows(dndb, is);
            // const auto dusdns0 = u0.ddn(is, is);

            // dns = dnsdb0 * (be - be0);
            // dus = dusdns0 * dns;
            // dus.noalias() = abs(dus / us0);

            uerror = dus.maxCoeff(&iuerror);

            if(uerror < reltol)
            {
                n = n0 + dndb * (be - be0);

                nmin = n.minCoeff(&inmin);
                ntot = sum(n);

                if(nmin/ntot < -1e-5)
                    continue;

                ranking[inode] += 1;

                auto comp = [&](Index l, Index r) { return ranking[l] > ranking[r]; };
                std::sort(priority.begin(), priority.end(), comp);

                // ChemicalState statetmp(state);
                // solver.solve(statetmp, T, P, be);

                // const auto utmp = solver.properties().chemicalPotentials().val;
                // const auto ntmp = statetmp.speciesAmounts();

                // Vector du = abs((utmp - u0.val)/u0.val);
                // Vector dn = abs((ntmp - n0)/n0);
                // auto actual_uerror = du.maxCoeff(&iuerror);
                // auto actual_nerror = dn.maxCoeff(&inerror);

                // if(actual_uerror > reltol && ntmp[iuerror] > 1e-10 && actual_nerror > reltol)
                // {
                //     std::cout << "==========================================================" << std::endl;
                //     std::cout << "species = " << system.species()[iuerror].name() << std::endl;
                //     std::cout << "uerror = " << uerror << std::endl;
                //     std::cout << "nerror = " << nerror << std::endl;
                //     std::cout << "n(exact) = " << ntmp[iuerror] << std::endl;
                //     std::cout << "n(trial) = " << n[iuerror] << std::endl;
                //     std::cout << "|n(exact) - n(trial)|/n(exact) = " << std::abs(n[iuerror] - ntmp[iuerror])/std::abs(ntmp[iuerror]) << std::endl;
                // }

                // const auto utmp = solver.properties().chemicalPotentials().val;
                // const auto ntmp = statetmp.speciesAmounts();

                // Vector du = abs((utmp - u0.val)/u0.val);
                // Vector dn = abs((ntmp - n0)/n0);
                // auto actual_uerror = du.maxCoeff(&iuerror);
                // auto actual_nerror = dn.maxCoeff(&inerror);

                // if(actual_uerror > reltol && ntmp[iuerror] > 1e-10 && actual_nerror > reltol)
                // {
                //     std::cout << "==========================================================" << std::endl;
                //     std::cout << "species = " << system.species()[iuerror].name() << std::endl;
                //     std::cout << "uerror = " << uerror << std::endl;
                //     std::cout << "nerror = " << nerror << std::endl;
                //     std::cout << "n(exact) = " << ntmp[iuerror] << std::endl;
                //     std::cout << "n(trial) = " << n[iuerror] << std::endl;
                //     std::cout << "|n(exact) - n(trial)|/n(exact) = " << std::abs(n[iuerror] - ntmp[iuerror])/std::abs(ntmp[iuerror]) << std::endl;
                // }

                // if(std::abs((n[n.size() - 1] - ntmp[n.size() - 1])/ntmp[n.size() - 1]) > 0.1 && ntmp[n.size() - 1] > 1e-10)
                // {
                //     std::cout << "--------------------------------------------------" << std::endl;
                //     std::cout << "Dolomite" << std::endl;
                //     std::cout << "n(exact) = " << ntmp[n.size() - 1] << std::endl;
                //     std::cout << "n(trial) = " << n[n.size() - 1] << std::endl;
                //     std::cout << "|n(exact) - n(trial)|/n(exact) = " << std::abs((n[n.size() - 1] - ntmp[n.size() - 1])/ntmp[n.size() - 1]) << std::endl;
                // }

                state.setSpeciesAmounts(n);
                // state.setElementDualPotentials(y);
                // state.setSpeciesDualPotentials(z);

                // Update the chemical properties of the system
                properties = node.properties;  // FIXME: We actually want to estimate properties = properties0 + variation : THIS IS A TEMPORARY SOLUTION!!!
                result.estimate.accepted = true;
                return;
            }
        }

        // std::cout << "ranking = " << std::endl;
        // for(auto i : priority)
        //     std::cout << ranking[i] << ", ";
        // std::cout << std::endl;

        if(result.estimate.accepted == false)
            return;







        MatrixConstRef Ae = partition.formulaMatrixEquilibriumPartition();

        // Comparison function based on the Euclidean distance
        auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        {
            Vector be_a = a.be/sum(a.be);
            Vector be_b = b.be/sum(b.be);
            Vector be_x = be/sum(be);

            return (be_a - be_x).squaredNorm() < (be_b - be_x).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        };

        // // Comparison function based on the Euclidean distance
        // auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        // {
        //     const auto& be_a = a.be;
        //     const auto& be_b = b.be;
        //     return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        // };

        // Comparison function based on the Euclidean distance
        // auto distancefn = [&](const TreeNode& a, const TreeNode& b)
        // {
        //     const auto RT = universalGasConstant*a.state.temperature();

        //     // Calculate perturbation of n
        //     a_dn.noalias() = a.sensitivity.dndb * (be - a.be);
        //     a_dy.noalias() = a.sensitivity.dydb * (be - a.be);
        //     a_dz.noalias() = a.sensitivity.dzdb * (be - a.be);

        //     a_n0.noalias() = a.state.speciesAmounts();
        //     a_y0.noalias() = a.state.elementDualPotentials();
        //     a_z0.noalias() = a.state.speciesDualPotentials();
        //     const auto& a_u0 = a.properties.chemicalPotentials();
        //     const auto& a_x0 = a.properties.moleFractions();

        //     a_n.noalias() = a_n0 + a_dn;
        //     a_y.noalias() = a_y0;
        //     a_z.noalias() = a_z0;


        //     b_dn.noalias() = b.sensitivity.dndb * (be - b.be);
        //     b_dy.noalias() = b.sensitivity.dydb * (be - b.be);
        //     b_dz.noalias() = b.sensitivity.dzdb * (be - b.be);

        //     b_n0.noalias() = b.state.speciesAmounts();
        //     b_y0.noalias() = b.state.elementDualPotentials();
        //     b_z0.noalias() = b.state.speciesDualPotentials();
        //     const auto& b_u0 = b.properties.chemicalPotentials();
        //     const auto& b_x0 = b.properties.moleFractions();

        //     b_n.noalias() = b_n0 + b_dn;
        //     b_y.noalias() = b_y0;
        //     b_z.noalias() = b_z0;


        //     // Correct negative mole numbers
        //     for(auto i = 0; i < n.size(); ++i)
        //     {
        //         // if(a_n[i] <= 0.0) a_n[i] = a_n0[i];
        //         // if(b_n[i] <= 0.0) b_n[i] = b_n0[i];
        //         if(a_n[i] <= 0.0) a_n[i] = 1e-12; // a_n0[i];
        //         if(b_n[i] <= 0.0) b_n[i] = 1e-12; // b_n0[i];
        //     }

        //     // Recompute the change in n
        //     a_dn = a_n - a_n0;
        //     b_dn = b_n - b_n0;

        //     // Calculate u(bar) = u(ref) + dudT(ref)*dT + dudP(ref)*dP + dudn(ref)*dn
        //     a_u.noalias() = a_u0.val + a_u0.ddn * a_dn;
        //     b_u.noalias() = b_u0.val + b_u0.ddn * b_dn;

        //     // Calculate x(bar) = x(ref) + dxdn(ref)*dn
        //     a_x.noalias() = a_x0.val + a_x0.ddn * a_dn;
        //     b_x.noalias() = b_x0.val + b_x0.ddn * b_dn;

        //     // Calculate the equilibrium residuals of the equilibrium species
        //     a_r.noalias() = abs(a_u - tr(Ae)*a_y - a_z)/RT;  // TODO: We should actually collect the entries in u and z corresponding to equilibrium species
        //     b_r.noalias() = abs(b_u - tr(Ae)*b_y - b_z)/RT;  // TODO: We should actually collect the entries in u and z corresponding to equilibrium species

        //     // // Eliminate species with mole fractions below cutoff from the residual analysis
        //     for(auto i = 0; i < n.size(); ++i)
        //     {
        //         if(a_x[i] < options.mole_fraction_cutoff) a_r[i] = 0.0; // set their residuals to zero
        //         if(b_x[i] < options.mole_fraction_cutoff) b_r[i] = 0.0; // set their residuals to zero
        //     }

        //     const double a_error = a_r.maxCoeff();
        //     const double b_error = b_r.maxCoeff();
        //     return a_error < b_error;

        //     // const auto& be_a = a.be;
        //     // const auto& be_b = b.be;
        //     // return (be_a - be).squaredNorm() < (be_b - be).squaredNorm();  // TODO: We need to extend this later with T and P contributions too (Allan, 26.06.2019)
        // };

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

        const auto RT = universalGasConstant*T0;

        // Calculate perturbation of n
        dn.noalias() = sensitivity0.dndb * (be - be0); // TODO: set derivatives dndb{i} = 0 when b{i} = 0
        dy.noalias() = sensitivity0.dydb * (be - be0); // TODO: set derivatives
        dz.noalias() = sensitivity0.dzdb * (be - be0); // TODO: set derivatives

        n.noalias() = n0 + dn;
        y.noalias() = y0;
        z.noalias() = z0;
        // y.noalias() = y0 + dy * RT; // TODO: Investigate further if derivatives of y and z wrt (T,P,b) can be made more accurately
        // z.noalias() = z0 + dz * RT; // TODO: Investigate further if derivatives of y and z wrt (T,P,b) can be made more accurately

        // Correct negative mole numbers
        for(auto i = 0; i < n.size(); ++i)
            if(n[i] <= 0.0)
                // n[i] = n0[i];
                n[i] = 1e-12;

        // Recompute the change in n
        dn = n - n0;

        // Calculate u(bar) = u(ref) + dudT(ref)*dT + dudP(ref)*dP + dudn(ref)*dn
        u = u0.val + u0.ddn * dn;

        // Calculate x(bar) = x(ref) + dxdn(ref)*dn
        x = x0.val + x0.ddn * dn;

        // Calculate the equilibrium residuals of the equilibrium species
        r = abs(u - tr(Ae)*y - z)/RT;  // TODO: We should actually collect the entries in u and z corresponding to equilibrium species

        // // Eliminate species with mole fractions below cutoff from the residual analysis
        for(auto i = 0; i < n.size(); ++i)
            if(x[i] < options.mole_fraction_cutoff)
                r[i] = 0.0; // set their residuals to zero

        // Eliminate species with amounts below their dual potentials (unstable species with n0[i] < z0[i])
        // for(auto i = 0; i < n.size(); ++i)
        //     if(n0[i] < z0[i])
        //         r[i] = 0.0; // set their residuals to zero





        // Should we add a test that prevent mole fractions from increasing by more than a certain tolerance?? My concern is that trace species could be left out.




        // Vector v = dn;

        // Eliminate species with mole fractions below cutoff from the residual analysis
        // for(auto i = 0; i < n.size(); ++i)
        //     if(x[i] < options.mole_fraction_cutoff)
        //         v[i] = 0.0; // set their residuals to zero

        // Eliminate species with amounts below their dual potentials (unstable species with n0[i] < z0[i])
        // for(auto i = 0; i < n.size(); ++i)
        //     if(n0[i] < z0[i])
        //         v[i] = 0.0; // set their residuals to zero

        // const double error = std::abs(v.dot(u0.ddn * v)/RT/sum(n0));
        // const double error = r.dot(abs(v))/sum(n);

        // Estimate the residual error of the trial Taylor approximation
        Index ispecies;
        const double error = r.maxCoeff(&ispecies);

        // const double error = abs(x.cwiseProduct(r)).maxCoeff(&ispecies);
        // const double error = r.maxCoeff(&ispecies); // find the maximum residual among considered stable and relevant species
        // abs(r.cwiseProduct(abs(v))/sum(n)).maxCoeff(&ispecies);
        // double error = 0.0;

        // Vector dx = abs(x - x0.val);

        // for(auto i = 0; i < n.size(); ++i)
        // {
        //     if(dx[i] > 1e-4) continue;
        //     if(dx[i]/x0.val[i] > 0.02) // 2% difference
        //     {
        //         ispecies = i;
        //         error = dx[i]/x0.val[i];
        //     }
        // }

        bool is_error_acceptable = error <= options.tol;
        // bool is_error_acceptable = true;

        if(is_error_acceptable == false)
        {
            std::cout << std::scientific;
            // std::cout << "-----------------------------" << std::endl;
            // std::cout << "*** Failed Smart Estimate ***" << std::endl;
            std::cout << "-----------------------------" << std::endl;
            std::cout << "Error = " << error << std::endl;
            std::cout << "Triggered by species = " << system.species(ispecies).name() << std::endl;
            std::cout << "-----------------------------" << std::endl;
            std::cout << std::left << std::setw(25) << "Species";
            std::cout << std::left << std::setw(25) << "r[i]";
            std::cout << std::left << std::setw(25) << "|r[i] * dn[i]/sum(n)|";
            std::cout << std::left << std::setw(25) << "|dn[i]|";
            std::cout << std::left << std::setw(25) << "n[i]";
            std::cout << std::left << std::setw(25) << "x[i]";
            std::cout << std::left << std::setw(25) << "z[i]";
            // std::cout << std::left << std::setw(25) << "r[i] * x[i]";
            // std::cout << std::left << std::setw(25) << "(r[i] * n[i])/nsum";
            std::cout << std::endl;
            for(auto i = 0; i < r.size(); ++i)
            {
                if(r[i] == 0.0) continue;
                if(i == ispecies) std::cout << "***************" << std::endl;
                std::cout << std::left << std::setw(25) << ((i == ispecies) ? "==> " : "") + system.species(i).name();
                std::cout << std::left << std::setw(25) << r[i];
                std::cout << std::left << std::setw(25) << r[i] * std::abs(dn[i])/sum(n);
                std::cout << std::left << std::setw(25) << std::abs(dn[i]);
                std::cout << std::left << std::setw(25) << n[i];
                std::cout << std::left << std::setw(25) << x[i];
                std::cout << std::left << std::setw(25) << z[i];
                // std::cout << std::left << std::setw(25) << std::abs(r[i] * x[i]);
                // std::cout << std::left << std::setw(25) << std::abs(r[i] * n[i])/nsum;
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

        state.setSpeciesAmounts(n);
        state.setElementDualPotentials(y);
        state.setSpeciesDualPotentials(z);

        toc(2, result.timing.estimate_acceptance);

        // Check if smart estimation failed with respect to variation of chemical potentials or amounts
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

