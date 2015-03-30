//// Reaktor is a C++ library for computational reaction modelling.
////
//// Copyright (C) 2014 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#include "KineticSolver.hpp"
//
//// C++ includes
//#include <functional>
//using namespace std::placeholders;
//
//// Eigen includes
//#include <eigen/LU>
//
//// Reaktor includes
//#include <Reaktor/Common/ChemicalVector.hpp>
//#include <Reaktor/Common/Constants.hpp>
//#include <Reaktor/Common/Exception.hpp>
//#include <Reaktor/Common/Matrix.hpp>
//#include <Reaktor/Common/Matrix.hxx>
//#include <Reaktor/Core/ChemicalState.hpp>
//#include <Reaktor/Core/ChemicalSystem.hpp>
//#include <Reaktor/Core/Partition.hpp>
//#include <Reaktor/Core/ReactionSystem.hpp>
//#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
//#include <Reaktor/Equilibrium/EquilibriumSolver.hpp>
//#include <Reaktor/Kinetics/KineticOptions.hpp>
//
//namespace Reaktor {
//
//class KineticSolver::Impl
//{
//private:
//    /// The chemical system instance
//    ChemicalSystem system;
//
//    /// The kinetically-controlled chemical reactions
//    ReactionSystem reactions;
//
//    /// The partition of the species in the chemical system
//    Partition partition;
//
//    /// The equilibrium solver instance
//    EquilibriumSolver equilibrium_solver;
//
//    /// The ODE solver instance
//    ODESolver ode_solver;
//
//    /// The ODE problem instance
//    ODEProblem ode_problem;
//
//    /// The evaluation of the right-hand side function of the ODE and its Jacobian
//    Vector f;
//
//    /// The evaluation of the Jacobian of the right-hand side function of the ODE
//    Matrix J;
//
//    /// The formula matrix of the equilibrium species
//    Matrix We;
//
//    /// The balance matrix of the equilibrium species
//    Matrix Be;
//
//    /// The kernel of the formula matrix of the equilibrium species
//    Matrix Ke;
//
//    /// The stoichiometric matrix w.r.t. the equilibrium species
//    Matrix Ve;
//
//    /// The stoichiometric matrix w.r.t. the kinetic species
//    Matrix Vk;
//
//    /// The elemental derivatives of the composition of the equilibrium species
//    Matrix Ze;
//
//    /// The molar abundance of the elements in the equilibrium species
//    Vector be;
//
//    /// The molar composition of the equilibrium species
//    Vector ne;
//
//    /// The molar composition of the kinetic species
//    Vector nk;
//
//    /// The combined vector of elemental molar abundance and composition of kinetic species [be nk]
//    Vector benk;
//
//    /// The activities of all species
//    ChemicalVector a;
//
//    /// The activities of the equilibrium species
//    ChemicalVector ae;
//
//    /// The kinetic rates of the reactions
//    ChemicalVector r;
//
//    /// The Jacobian of the kinetic rate w.r.t. the equilibrium species
//    Matrix Re;
//
//    /// The Jacobian of the kinetic rate w.r.t. the kinetic species
//    Matrix Rk;
//
//    /// The coefficient matrix of the equations of the elemental derivatives of the equilibrium species
//    Matrix He;
//
//    /// The partial pivoting LU instance for the solution of the elemental derivatives of the equilibrium species
//    Eigen::PartialPivLU<Matrix> lu;
//
//    /// The number of equilibrium species in the system
//    unsigned Ne;
//
//    /// The number of kinetic species in the system
//    unsigned Nk;
//
//    /// The number of elements in the system
//    unsigned Nb;
//
//    /// The number of constraints in the system (mass-balance and charge-balance equations)
//    unsigned Nc;
//
//    /// The universal gas constant times the temperature in the system
//    double RT;
//
//public:
//    Impl(const ChemicalSystem& system, const Partition& partition, const ReactionSystem& reactions)
//    : system(system), reactions(reactions), partition(partition),
//      equilibrium_solver(system, partition)
//    {
//        BalanceConstraints balance(system, partition);
//
//        // Set the number of equilibrium and kinetic species, and the number of elements
//        Ne = partition.numEquilibriumSpecies();
//        Nk = partition.numKineticSpecies();
//        Nb = partition.numEquilibriumElements();
//        Nc = balance.numConstraints();
//
//        // Initialise the formula matrix of the equilibrium species
//        We = partition.equilibriumFormulaMatrix(system);
//
//        // Initialise the balance matrix of the equilibrium species
//        Be = balance.balanceMatrix();
//
//        // Initialise the kernel of the formula matrix of the equilibrium species
//        Ke = Be.fullPivLu().kernel();
//
//        // Initialise the elemental derivatives of the composition of the equilibrium species
//        Ze = zeros(Ne, Nb);
//
//        // Initialise the stoichiometric matrices w.r.t. the equilibrium and kinetic species
//        Ve = partition.equilibriumStoichiometricMatrix(reactions);
//        Vk = partition.kineticStoichiometricMatrix(reactions);
//
//        // Initialise the result of the evaluation of the ODE function and its Jacobian
//        f.resize(Nb + Nk);
//        J.resize(Nb + Nk, Nb + Nk);
//
//        // Set the ODE problem
//        ode_problem.setNumEquations(Nb + Nk);
//        ode_problem.setFunction(std::bind(&Impl::function, this, _1, _2, _3));
//        ode_problem.setJacobian(std::bind(&Impl::jacobian, this, _1, _2, _3));
//
//        // Set the ODE problem in the ODE solver
//        ode_solver.setProblem(ode_problem);
//    }
//
//    auto setOptions(const KineticOptions& options) -> void
//    {
//        ode_solver.setOptions(options.ode);
//        equilibrium_solver.setOptions(options.equilibrium);
//    }
//
//    auto setProblem(const KineticProblem& _problem) -> void
//    {
//        problem = _problem;
//    }
//
//    auto initialise(const ChemicalState& state, double tstart) -> void
//    {
//        // The temperature of the chemical system
//        const double T = state.temperature();
//
//        // Extract the composition vector of the equilibrium and kinetic species
//        ne = state.equilibriumComposition(partition);
//        nk = state.kineticComposition(partition);
//
//        // Set the universal gas constant times the temperature variable
//        RT = universalGasConstant * T;
//
//        // Assemble the vector benk = [be nk]
//        benk = zeros(Nb + Nk);
//        benk.segment(00, Nb) = We * ne;
//        benk.segment(Nb, Nk) = nk;
//
//        // Initialise the ODE solver
//        ode_solver.Initialise(tstart, benk);
//    }
//
//    auto step(ChemicalState& state, double& t) -> void
//    {
//        const double tfinal = unsigned(-1);
//        step(state, t, tfinal);
//    }
//
//    auto step(ChemicalState& state, double& t, double tfinal) -> void
//    {
//        // Set the state member of the equilibrium_result instance to be used by Function and Jacobian
//        state = state;
//
//        // Extract the composition vector of the equilibrium and kinetic species
//        ne = state.equilibriumComposition(partition);
//        nk = state.kineticComposition(partition);
//
//        // Assemble the vector benk = [be nk]
//        benk.segment(00, Nb) = We * ne;
//        benk.segment(Nb, Nk) = nk;
//
//        // Perform one ODE step integration
//        ode_solver.Integrate(t, benk, tfinal);
//
//        // Extract the `be` and `nk` entries of the vector `benk`
//        be = benk.segment(00, Nb);
//        nk = benk.segment(Nb, Nk);
//
//        // Update the composition of the kinetic species
//        state.setKineticComposition(partition, nk);
//
//        // Update the composition of the equilibrium species
//        equilibrium_solver.minimise(state, be);
//
//        // Update the provided chemical state with the internal one
//        state = state;
//
//        lagrange = lagrange;
//    }
//
//    auto solve(ChemicalState& state, double tstart, double tfinal) -> void
//    {
//        // Extract the composition vector of the equilibrium and kinetic species
//        ne = state.equilibriumComposition(partition);
//        nk = state.kineticComposition(partition);
//
//        // Assemble the vector benk = [be nk]
//        benk = zeros(Nb + Nk);
//        benk.segment(00, Nb) = We * ne;
//        benk.segment(Nb, Nk) = nk;
//
//        // Perform one ODE step integration
//        ode_solver.Solve(tstart, tfinal, benk);
//
//        // Extract the `be` and `nk` entries of the vector `benk`
//        be = benk.segment(00, Nb);
//        nk = benk.segment(Nb, Nk);
//
//        // Update the composition of the kinetic species
//        state.setKineticComposition(partition, nk);
//
//        // Update the composition of the equilibrium species
//        equilibrium_solver.minimise(state, lagrange, be);
//    }
//
//    auto function(double t, const Vector& benk, Vector& f) -> int
//    {
//        // Extract the `be` and `nk` entries of the vector [be, nk]
//        be = benk.segment(00, Nb);
//        nk = benk.segment(Nb, Nk);
//
//        // Check for non-finite values in the vector `benk`
//        for(unsigned i = 0; i < benk.rows(); ++i)
//            if(not std::isfinite(benk[i]))
//                return 1; // ensure the ode solver will reduce the time step
//
//        // Do not allow that the molar amount of the elements be negative
//        if(be.minCoeff() < 0.0)
//            return 1; // ensure the ode solver will reduce the time step
//
//        // Update the composition of the kinetic species in the member `state`
//        state.setKineticComposition(partition, nk);
//
//        // Solve the equilibrium problem using the elemental molar abundance `be`
//        equilibrium_solver.minimise(state, lagrange, be);
//
//        // Use the just calculated system composition to update the activities of the species
//        a = state.activities();
//
//        // Calculate the kinetic rates of the reactions
//        r = reactions.rates(state, a);
//
//        // Calculate the right-hand side function of the ODE
//        f.segment(00, Nb) = We * Ve.transpose() * func(r);
//        f.segment(Nb, Nk) = Vk.transpose() * func(r);
//
//        // Impose a lower bound for the decrease of some kinetic species
//        for(unsigned i = 0; i < benk.rows(); ++i)
//            if(std::abs(benk[i]) < 1.0e-50 and f[i] < 0.0)
//                f[i] = 0.0; // set the rate to zero
//
//        return 0;
//    }
//
//    auto jacobian(double t, const Vector& benk, Matrix& J) -> int
//    {
//        // Extract the `be` and `nk` entries of the vector [be, nk]
//        be = benk.segment(00, Nb);
//        nk = benk.segment(Nb, Nk);
//
//        // Update the composition of the kinetic species in the member `state`
//        state.setKineticComposition(partition, nk);
//
//        // Solve the equilibrium problem using the elemental molar abundance `be`
//        equilibrium_solver.minimise(state, lagrange, be);
//
//        // Use the just calculated system composition to update the activities of the species
//        a = state.activities();
//
//        // Extract the activities and its derivatives of the equilibrium species
//        func(ae) = partition.equilibriumRows(func(a));
//        grad(ae) = partition.equilibriumRowsCols(grad(a));
//
//        // The Lagrange multipliers `ze` of the equilibrium calculation
//        Vector ze = lagrange.lagrangeZ();
//
//        // Assemble the coefficient matrix of the elemental derivative equations
//        He.resize(Ne, Ne);
//        He.topRows(Ne - Nc)  = Ke.transpose() * (func(ae).asDiagonal().inverse() * grad(ae));
//        He.topRows(Ne - Nc) += Ke.transpose() * ze.cwiseQuotient(ne).asDiagonal();
//        He.bottomRows(Nc)    = Be;
//
//        // Perfom a LU decomposition of the coefficient matrix
//        lu.compute(He);
//
//        // Check for non-finite values in the LU decomposition
//        if(lu.matrixLU() != lu.matrixLU())
//            RuntimeError("Cannot compute the jacobian matrix.", "There is a nan in the LU(He).");
//
//        // Calculate the right-hand side matrix
//        Matrix rhs = zeros(Ne, Nb);
//        rhs.middleRows(Ne - Nc, Nb).setIdentity();
//
//        // Calculate the elemental derivatives of the composition of the equilibrium species
//        Ze = lu.solve(rhs);
//
//        // Extract the equilibrium and kinetic entries of the Jacobian of the kinetic rates
//        Re = partition.equilibriumCols(grad(r));
//        Rk = partition.kineticCols(grad(r));
//
//        // Calculate the Jacobian of the right-hand side function of the ODE
//        J.block(00, 00, Nb, Nb) = We * Ve.transpose() * Re * Ze;
//        J.block(00, Nb, Nb, Nk) = We * Ve.transpose() * Rk;
//        J.block(Nb, 00, Nk, Nb) = Vk.transpose() * Re * Ze;
//        J.block(Nb, Nb, Nk, Nk) = Vk.transpose() * Rk;
//
//        return 0;
//    }
//};
//
//KineticSolver::KineticSolver()
//: pimpl(new Impl())
//{}
//
//KineticSolver::KineticSolver(const ReactionSystem& reactions)
//: pimpl(new Impl(reactions))
//{}
//
//KineticSolver::KineticSolver(const KineticSolver& other)
//: pimpl(new Impl(*other.pimpl))
//{}
//
//KineticSolver::~KineticSolver()
//{}
//
//auto KineticSolver::operator=(KineticSolver other) -> KineticSolver&
//{
//    pimpl = std::move(other.pimpl);
//    return *this;
//}
//
//auto KineticSolver::setOptions(const KineticOptions& options) -> void
//{
//    pimpl->setOptions(options);
//}
//
//auto KineticSolver::setProblem(const KineticProblem& problem) -> void
//{
//    pimpl->setProblem(problem);
//}
//
//auto KineticSolver::step(ChemicalState& state, double& t) -> void
//{
//    pimpl->step(state, t);
//}
//
//auto KineticSolver::step(ChemicalState& state, double& t, double dt) -> void
//{
//    pimpl->step(state, t, dt);
//}
//
//auto KineticSolver::solve(ChemicalState& state, double t, double dt) -> void
//{
//    pimpl->solve(state, tstart, tfinal);
//}
//
//} // namespace Reaktor
