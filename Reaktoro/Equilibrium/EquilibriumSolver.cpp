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

#include "EquilibriumSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Connectivity.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Math/MathUtils.hpp>
#include <Reaktoro/Optimization/OptimumOptions.hpp>
#include <Reaktoro/Optimization/OptimumProblem.hpp>
#include <Reaktoro/Optimization/OptimumResult.hpp>
#include <Reaktoro/Optimization/OptimumSolver.hpp>
#include <Reaktoro/Optimization/OptimumSolverRefiner.hpp>
#include <Reaktoro/Optimization/OptimumState.hpp>

#include <Reaktoro/Thermodynamics/Core/Database.hpp>

phaseIdentificationMethod phaseidMethod;
extern int quantity;

namespace Reaktoro {

struct EquilibriumSolver::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The options of the equilibrium solver
    EquilibriumOptions options;

    /// The solver for the optimisation calculations
    OptimumSolver solver;

    /// The chemical properties of the chemical system
    ChemicalProperties properties;

    /// The sensitivity derivatives of the equilibrium state
    EquilibriumSensitivity sensitivities;
    Vector zerosEe; // FIXME: Improve design. These vectors are needed to calculate sensitivities, but they should not exist!
    Vector zerosNe; // FIXME: Improve design. These vectors are needed to calculate sensitivities, but they should not exist!
    Vector unitjEe; // FIXME: Improve design. These vectors are needed to calculate sensitivities, but they should not exist!

    /// The molar amounts of the species
    Vector n;

    /// The dual potentials of the elements
    Vector y;

    /// The dual potentials of the species
    Vector z;

    /// The molar amounts of the elements in the equilibrium partition
    Vector be;

    /// The standard chemical potentials of the species
    ThermoVector u0;

    /// The chemical potentials of the species
    ChemicalVector u;

    /// The chemical potentials of the equilibrium species
    ChemicalVector ue;

    /// The chemical potentials of the inert species
    Vector ui;

    /// The mole fractions of the equilibrium species
    ChemicalVector xe;

    /// The optimisation problem
    OptimumProblem optimum_problem;

    /// The state of the optimisation calculation
    OptimumState optimum_state;

    /// The options for the optimisation calculation
    OptimumOptions optimum_options;

    /// The indices of the species in the equilibrium partition
    Indices ies;

    /// The indices of the elements in the equilibrium partition
    Indices iee;

    /// The indices of the inert species (i.e., the species in disequilibrium)
    Indices iis;

    /// The number of species and elements in the system
    unsigned N, E;

    /// The number of species and elements in the equilibrium partition
    unsigned Ne, Ee;

    /// The formula matrix of the species in the system
    Matrix A;

    /// The formula matrix of the species in the equilibrium partition
    Matrix Ae;

    /// The formula matrix of the inert species
    Matrix Ai;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance
    Impl(const ChemicalSystem& system)
    : system(system), properties(system)
    {
        // Initialize the formula matrix
        A = system.formulaMatrix();

        // Initialize the number of species and elements in the system
        N = system.numSpecies();
		E = system.numElements();

		// Set the default partition as all species are in equilibrium
		setPartition(Partition(system));
	}

	/// Set the partition of the chemical system
	auto setPartition(const Partition& partition_) -> void
	{
		// Set the partition of the chemical system
		partition = partition_;

		// Initialize the number of species and elements in the equilibrium partition
		Ne = partition.numEquilibriumSpecies();
		Ee = partition.numEquilibriumElements();

		// Initialize the formula matrix of the equilibrium species
		Ae = partition.formulaMatrixEquilibriumPartition();

		// Initialize the indices of the equilibrium species and elements
		ies = partition.indicesEquilibriumSpecies();
		iee = partition.indicesEquilibriumElements();

		// Initialize the indices of the inert species
		iis.clear();
		iis.reserve(partition.numInertSpecies() + partition.numKineticSpecies());
		iis.insert(iis.end(), partition.indicesInertSpecies().begin(), partition.indicesInertSpecies().end());
		iis.insert(iis.end(), partition.indicesKineticSpecies().begin(), partition.indicesKineticSpecies().end());

		// Initialize the formula matrix of the inert species
		Ai = cols(A, iis);
	}

	/// Update the OptimumOptions instance with given EquilibriumOptions instance
	auto updateOptimumOptions() -> void
	{
		// Initialize the options for the optimisation calculation
		optimum_options = options.optimum;

		// Set the parameters of the optimisation algorithms that control how small can be the amount of a species
		optimum_options.ipaction.mu = options.epsilon;
		optimum_options.ipnewton.mu = options.epsilon;
		optimum_options.ipopt.mu.push_back(options.epsilon);
		optimum_options.ipactive.epsilon = options.epsilon;

		// Initialize the names of the primal and dual variables
		if (options.optimum.output.active)
		{
			// Use `n` instead of `x` to name the variables
			optimum_options.output.xprefix = "n";

			// Define some auxiliary references to the variables names
			auto& xnames = optimum_options.output.xnames;
			auto& ynames = optimum_options.output.ynames;
			auto& znames = optimum_options.output.znames;

			// Initialize the names of the primal variables `n`
			for (Index i : ies)
				xnames.push_back(system.species(i).name());

			// Initialize the names of the dual variables `y`
			for (Index i : iee)
				ynames.push_back(system.element(i).name());

			// Initialize the names of the dual variables `z`
			znames = xnames;
		}
	}

	/// Update the OptimumProblem instance with given EquilibriumProblem and ChemicalState instances
	auto updateOptimumProblem(const ChemicalState& state) -> void
	{
		// The temperature and pressure of the equilibrium calculation
		const auto T = state.temperature();
		const auto P = state.pressure();
		const auto RT = universalGasConstant * T;

		// Set the molar amounts of the species
		n = state.speciesAmounts();

		// Update the standard thermodynamic properties of the chemical system
		properties.update(T, P);

		// Update the normalized standard Gibbs energies of the species
		u0 = properties.standardPartialMolarGibbsEnergies() / RT;

		// The result of the objective evaluation
		ObjectiveResult res;

		// The Gibbs energy function to be minimized
		optimum_problem.objective = [=](VectorConstRef ne) mutable
		{
			// Set the molar amounts of the species
			n(ies) = ne;

			// Update the chemical properties of the chemical system
			properties.update(T, P, n);

            // Set the scaled chemical potentials of the species
            u = u0 + properties.lnActivities();

            // Set the scaled chemical potentials of the equilibrium species
            ue = rows(u, ies, ies);

            // Set the mole fractions of the equilibrium species
            xe = rows(properties.moleFractions(), ies, ies);

            // Set the objective result
            res.val = dot(ne, ue.val);
            res.grad = ue.val;

            // Set the Hessian of the objective function
            switch(options.hessian)
            {
            case GibbsHessian::Exact:
                res.hessian.mode = Hessian::Dense;
                res.hessian.dense = ue.ddn;
                break;
            case GibbsHessian::ExactDiagonal:
                res.hessian.mode = Hessian::Diagonal;
                res.hessian.diagonal = diagonal(ue.ddn);
                break;
            case GibbsHessian::Approximation:
                res.hessian.mode = Hessian::Dense;
                res.hessian.dense = diag(inv(xe.val)) * xe.ddn;
                break;
            case GibbsHessian::ApproximationDiagonal:
                res.hessian.mode = Hessian::Diagonal;
                res.hessian.diagonal = diagonal(xe.ddn)/xe.val;
                break;
            }

            return res;
        };

        optimum_problem.c.resize(0);
        optimum_problem.n = Ne;
        optimum_problem.A = Ae;
        optimum_problem.b = be;
        optimum_problem.l.setConstant(Ne, options.epsilon);
    }

    /// Initialize the optimum state from a chemical state
    auto updateOptimumState(const ChemicalState& state) -> void
    {
        // The temperature and the RT factor
        const double T  = state.temperature();
        const double RT = universalGasConstant*T;

        // Set the molar amounts of the species
        n = state.speciesAmounts();

        // Set the normalized dual potentials of the elements
        y = state.elementDualPotentials()/RT;

        // Set the normalized dual potentials of the species
        z = state.speciesDualPotentials()/RT;

        // Initialize the optimum state
        optimum_state.x = n(ies);
        optimum_state.y = y(iee);
        optimum_state.z = z(ies);
    }

    /// Initialize the chemical state from a optimum state
    auto updateChemicalState(ChemicalState& state) -> void
    {
        // The temperature and the RT factor
        const double T  = state.temperature();
        const double RT = universalGasConstant*T;

        // Update the molar amounts of the equilibrium species
        n(ies) = optimum_state.x;

        // Update the normalized chemical potentials of the inert species
        ui = u.val(iis);

        // Update the normalized dual potentials of the elements
        y = zeros(E); y(iee) = optimum_state.y;

        // Update the normalized dual potentials of the equilibrium and inert species
        z(ies) = optimum_state.z;
        z(iis) = ui - tr(Ai) * y;

        // Scale the normalized dual potentials of elements and species to units of J/mol
        y *= RT;
        z *= RT;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setElementDualPotentials(y);
        state.setSpeciesDualPotentials(z);
    }

    /// Find a feasible approximation for an equilibrium problem.
    auto approximate(ChemicalState& state, double T, double P, Vector be) -> EquilibriumResult
    {
        // Check the dimension of the vector `be`
        Assert(unsigned(be.rows()) == Ee,
            "Cannot proceed with method EquilibriumSolver::approximate.",
            "The dimension of the given vector of molar amounts of the "
            "elements does not match the number of elements in the "
            "equilibrium partition.");

        // Set temperature and pressure of the chemical state
        state.setTemperature(T);
        state.setPressure(P);

        // Auxiliary variables
        const double RT = universalGasConstant*T;
        const double inf = std::numeric_limits<double>::infinity();

        // Update the internal state of n, y, z
        n = state.speciesAmounts();
        y = state.elementDualPotentials();
        z = state.speciesDualPotentials();

        // Update the standard thermodynamic properties of the system
        properties.update(T, P);

        // Get the standard Gibbs energies of the equilibrium species
        const Vector ge0 = properties.standardPartialMolarGibbsEnergies().val(ies);

        // Get the ln activity constants of the equilibrium species
        const Vector ln_ce = properties.lnActivityConstants().val(ies);

        // Define the optimisation problem
        OptimumProblem optimum_problem;
        optimum_problem.n = Ne;
        optimum_problem.c = ge0/RT + ln_ce;
        optimum_problem.A = Ae;
        optimum_problem.b = be;
        optimum_problem.l = zeros(Ne);
        optimum_problem.u = ones(Ne) * inf;

        // Initialize the optimum state
        OptimumState optimum_state;

        // The result of the linear programming calculation
        EquilibriumResult result;

        // Update the optimum options
        updateOptimumOptions();

        // Set the method for the optimisation calculation
        solver.setMethod(options.method);

        // Solve the linear programming problem
        result.optimum = solver.solve(optimum_problem, optimum_state, optimum_options);

        // Update the molar amounts of the equilibrium species
        n(ies) = optimum_state.x;

        // Update the dual potentials of the species and elements (in units of J/mol)
        z = zeros(N); z(ies) = optimum_state.z * RT;
        y = zeros(E); y(iee) = optimum_state.y * RT;

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setElementDualPotentials(y);
        state.setSpeciesDualPotentials(z);

        return result;
    }

    /// Find an initial guess for an equilibrium problem.
    auto initialguess(ChemicalState& state, double T, double P, Vector be) -> EquilibriumResult
    {
        // Solve the linear programming problem to obtain an approximation
        auto result = approximate(state, T, P, be);

        // Check the approximate calculation was successful
        Assert(result.optimum.succeeded,
            "Cannot proceed with the equilibrium calculation.",
            "The calculation of initial guess failed, most "
            "probably because no feasible solution exists.");

        // Preserve the values of n that are greater than z. For all others, set n to sqrt(epsilon)
        n = (n.array() > z.array()).select(n, std::sqrt(options.epsilon));

        // Set z to zero, which will later be set to epsilon / n.
        z.fill(0.0);

        // Update the chemical state
        state.setSpeciesAmounts(n);
        state.setSpeciesDualPotentials(z);

        return result;
    }

    /// Return true if cold-start is needed.
    auto coldstart(const ChemicalState& state) -> bool
    {
        // Check if all equilibrium species have zero amounts
        bool zero = true;
        for(Index i : ies)
            if(state.speciesAmount(i) > 0)
                { zero = false; break; }
        return zero || !options.warmstart;
    }
    /*
	/// Return true if the system has the need to run the phase identification to remove an inaproriet phase
	auto necessityofphaseidentification(const ChemicalState& state) -> bool
	{
		// Check the presence of Gas and Oil
		const auto& phases = state.system().phases();
		auto num_oil_and_gas = 0;
		auto number_possible_phases = phases.size();
		auto z = Vector(4);
		for (const auto& phase : phases) {
			if (phase.name() == "Oil" || phase.name() == "Gaseous"){
				num_oil_and_gas++;
			}
		}
		if (num_oil_and_gas > 1) {
			z = speciesAmountsNotaqueous(state);
			auto system = ChemicalSystem({ state.system().phase("Gaseous") });
			number_possible_phases = stabitytest_TPD(system, z, state.temperature(), state.pressure());
			if (number_possible_phases <= 1) {
				return true;
			}
		}
		return false;
	}
    */
	
	/// Return a vector with the total amount of species without count Aqueous phase
	// the currently version only works when all species in oil phase are also in gas
	auto speciesAmountsNotaqueous(const ChemicalState& state) -> Vector
	{
		const auto& system = state.system();
		auto z = Vector(4);
		const auto gas_phase_species = system.phase("Gaseous").species();
		for (auto i = 0; i < gas_phase_species.size(); i++) {
			if (gas_phase_species[i].name() == "H2O(g)") {
				z[i] = state.speciesAmount(gas_phase_species[i].name());
			}
			else {
				auto found = (gas_phase_species[i].name()).find("(g)");
				auto gas_specie_name = gas_phase_species[i].name().substr(0, found) + "(g)";
				auto oil_specie_name = gas_phase_species[i].name().substr(0, found) + "(oil)";
				z[i] = state.speciesAmount(gas_specie_name) + state.speciesAmount(oil_specie_name);
			}
		}
		return z;
	}

    /*
	/// Return the number os possible stable phase of the system and the vector of element amounts
	auto stabitytest_TPD(const ChemicalSystem& system, Vector Z, double temperature, double pressure) -> int {
		auto numberofspecies = system.species().size();
		auto K_liq = Vector(numberofspecies);
		auto K_vap = Vector(numberofspecies);
		auto Y_liq = Vector(numberofspecies);
		auto y_liq = Vector(numberofspecies);
		auto Y_vap = Vector(numberofspecies);
		auto y_vap = Vector(numberofspecies);

		auto R_liq = Vector(numberofspecies);
		auto R_vap = Vector(numberofspecies);

		auto db = Database("W:\\release\\Projects\\Reaktoro\\databases\\supcrt\\supcrt98.xml");

		auto z_cres = ChemicalModelResult(numberofspecies, numberofspecies);
		
		for (auto i = 0; i < Z.size(); i++) {
			z[i] = z[i] / Z.sum();
		}
		
		system.chemicalModel()(z_cres, temperature, pressure, z);


		for (auto i = 0; i < system.numSpecies(); i++) {							
			double Pc = db.gaseousSpecies(system.species(i).name()).criticalPressure();
			double Tc = db.gaseousSpecies(system.species(i).name()).criticalTemperature();
			double w = db.gaseousSpecies(system.species(i).name()).acentricFactor();
			auto Ki = (Pc / pressure)*std::exp(5.37 * (1.0 + w)*(1.0 - (Tc / temperature)));
			K_liq(i) = Ki;
			K_vap(i) = Ki;
		}

		bool converged = false;
		bool TSliq = false;
		bool TSvap = false;

		double Sl = 0.0;
		double Sv = 0.0;
		for (auto i = 0; i < 10000; i++) {
			for (auto j = 0; j < system.numSpecies(); j++) {
				Y_liq(j) = z(j) / K_liq(j);
				Y_vap(j) = z(j)*K_vap(j);
			}

			Sl = Y_liq.sum();
			Sv = Y_vap.sum();

			for (auto j = 0; j < system.numSpecies(); j++) {
				y_liq(j) = Y_liq(j) / Sl;
				y_vap(j) = Y_vap(j) / Sv;
			}

			auto y_liq_cres = ChemicalModelResult(numberofspecies, numberofspecies);
			auto y_vap_cres = ChemicalModelResult(numberofspecies, numberofspecies);

			system.chemicalModel()(y_liq_cres, temperature, pressure, y_vap);
			system.chemicalModel()(y_vap_cres, temperature, pressure, y_liq);

			for (auto j = 0; j < system.numSpecies(); j++) {
				R_liq(j) = (Sl)*(std::exp(y_liq_cres.lnActivityCoefficients()[j].val) / std::exp(z_cres.lnActivityCoefficients()[j].val));
				R_vap(j) = (1.0 / Sv)*(std::exp(z_cres.lnActivityCoefficients()[j].val) / std::exp(y_vap_cres.lnActivityCoefficients()[j].val));
			}

			auto trial_solution_liq = 0.0;
			auto trial_solution_vap = 0.0;
			for (auto j = 0; j < system.numSpecies(); j++) {
				trial_solution_liq += std::log(K_liq(j))*std::log(K_liq(j));
				trial_solution_vap += std::log(K_vap(j))*std::log(K_vap(j));
			}
			if (trial_solution_liq < 1.e-4) {
				TSliq = true;
			}
			if (trial_solution_vap < 1.e-4) {
				TSvap = true;
			}

			auto erro_liq = 0.0;
			auto erro_vap = 0.0;
			for (auto j = 0; j < system.numSpecies(); j++) {
				erro_liq += (R_liq(j) - 1)*(R_liq(j) - 1);
				erro_vap += (R_vap(j) - 1)*(R_vap(j) - 1);
			}
			if (erro_liq < 1e-12) {
				//std::cout << "liq possivel to form" << std::endl;
				converged = true;
			}
			if (erro_vap < 1e-12) {
				//std::cout << "vap to form" << std::endl;
				converged = true;
			}
			if (converged) {
				break;
			}

			for (auto j = 0; j < system.numSpecies(); j++) {
				K_liq(j) = K_liq(j) * R_liq(j);
				K_vap(j) = K_vap(j) * R_liq(j);
			}
			if (i == 99999) {
				return 4;
			}
		}

		if (converged) {
			if (TSliq && TSvap) {
				return 1; //std::cout << "1 fase" << std::endl;
			}
			else if (Sl <= 1.0 && TSvap) {
				return 2; // std::cout << "2 fases" << std::endl;
			}
			else if (TSliq && Sl <= 1.0) {
				return 2;// std::cout << "2 fases" << std::endl;
			}
			else if (Sv <= 1.0 && Sl <= 1.0) {
				return 3;// std::cout << "3 fases" << std::endl;
			}
		}
		else {
			if (Sv > 1 && TSliq) {
				return 1;// std::cout << "2 fases" << std::endl;
			}
			else if (TSvap && Sl > 1) {
				return 2;// std::cout << "2 fases" << std::endl;
			}
			else if (Sv > 1 && Sl > 1) {
				return 2;// std::cout << "2 fases" << std::endl;
			}
			else if (Sv > 1 && Sl <= 1) {
				return 3;// std::cout << "3 fases" << std::endl;
			}
			else if (Sv <= 1 && Sl > 1) {
				return 3;// std::cout << "3 fases" << std::endl;
			}
		}
		return 4;
	}
    */
    /// Solve the equilibrium problem
    auto solve(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
    {
        // Check the dimension of the vector `be`
        Assert(be.size() == static_cast<int>(Ee),
            "Cannot proceed with method EquilibriumSolver::solve.",
            "The dimension of the given vector of molar amounts of the "
            "elements does not match the number of elements in the "
            "equilibrium partition.");
        return solve(state, T, P, be.data());
    }

    /// Solve the equilibrium problem
    auto solve(ChemicalState& state, double T, double P, const double* b) -> EquilibriumResult
    {
        // Set the molar amounts of the elements
        be = Vector::Map(b, Ee);

        // Set temperature and pressure of the chemical state
        state.setTemperature(T);
        state.setPressure(P);

        // Check if a simplex cold-start approximation must be performed
        if(coldstart(state))
            initialguess(state, T, P, be);

        // The result of the equilibrium calculation
        EquilibriumResult result;

        // Update the optimum options
        updateOptimumOptions();

        // Update the optimum problem
        updateOptimumProblem(state);

        // Update the optimum state
        updateOptimumState(state);

        // Set the method for the optimisation calculation
        solver.setMethod(options.method);

        // Set the maximum number of iterations in each optimization pass
        optimum_options.max_iterations = 10;

        // Start the several opmization passes (stop if convergence attained)
        auto counter = 0;
        while(counter < options.optimum.max_iterations)
        {
            // Solve the optimisation problem
            result.optimum += solver.solve(optimum_problem, optimum_state, optimum_options);

            // Exit this loop if last solve succeeded
            if(result.optimum.succeeded)
                break;

            counter += optimum_options.max_iterations;
        }

        // Update the chemical state from the optimum state
        updateChemicalState(state);

        return result;
    }

    /// Return the sensitivity of the equilibrium state.
    auto sensitivity() -> const EquilibriumSensitivity&
    {
        zerosEe = zeros(Ee);
        zerosNe = zeros(Ne);
        unitjEe = zeros(Ee);

        sensitivities.dndT = zeros(Ne);
        sensitivities.dndP = zeros(Ne);
        sensitivities.dndb = zeros(Ne, Ee);

        sensitivities.dndT = solver.dxdp(ue.ddT, zerosEe);
        sensitivities.dndP = solver.dxdp(ue.ddP, zerosEe);
        for(Index j = 0; j < Ee; ++j)
        {
            unitjEe = unit(Ee, j);
            sensitivities.dndb.col(j) = solver.dxdp(zerosNe, unitjEe);
        }

        return sensitivities;
    }

    /// Compute the sensitivity of the species amounts with respect to temperature.
    auto dndT() -> VectorConstRef
    {
        const auto& ieq_species = partition.indicesEquilibriumSpecies();
        zerosEe = zeros(Ee);
        sensitivities.dndT = zeros(N);
        sensitivities.dndT(ieq_species) = solver.dxdp(ue.ddT, zerosEe);
        return sensitivities.dndT;
    }

    /// Compute the sensitivity of the species amounts with respect to pressure.
    auto dndP() -> VectorConstRef
    {
        const auto& ieq_species = partition.indicesEquilibriumSpecies();
        zerosEe = zeros(Ee);
        sensitivities.dndP = zeros(N);
        sensitivities.dndP(ieq_species) = solver.dxdp(ue.ddP, zerosEe);
        return sensitivities.dndP;
    }

    /// Compute the sensitivity of the species amounts with respect to element amounts.
    auto dndb() -> VectorConstRef
    {
        const auto& ieq_species = partition.indicesEquilibriumSpecies();
        const auto& ieq_elements = partition.indicesEquilibriumElements();
        zerosEe = zeros(Ee);
        zerosNe = zeros(Ne);
        unitjEe = zeros(Ee);
        sensitivities.dndb = zeros(Ne, Ee);
        for(Index j : ieq_elements)
        {
            unitjEe = unit(Ee, j);
            sensitivities.dndb.col(j)(ieq_species) = solver.dxdp(zerosNe, unitjEe);
        }
        return sensitivities.dndb;
    }
};

EquilibriumSolver::EquilibriumSolver()
: pimpl(new Impl())
{}

EquilibriumSolver::EquilibriumSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

EquilibriumSolver::EquilibriumSolver(const EquilibriumSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

EquilibriumSolver::~EquilibriumSolver()
{}

auto EquilibriumSolver::operator=(EquilibriumSolver other) -> EquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto EquilibriumSolver::setOptions(const EquilibriumOptions& options) -> void
{
    pimpl->options = options;
}

auto EquilibriumSolver::setPartition(const Partition& partition) -> void
{
    pimpl->setPartition(partition);
}

auto EquilibriumSolver::approximate(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
{
    return pimpl->approximate(state, T, P, be);
}

auto EquilibriumSolver::approximate(ChemicalState& state, const EquilibriumProblem& problem) -> EquilibriumResult
{
    return approximate(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto EquilibriumSolver::approximate(ChemicalState& state) -> EquilibriumResult
{
    return approximate(state, state.temperature(), state.pressure(), state.elementAmounts());
}

auto EquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto EquilibriumSolver::solve(ChemicalState& state, double T, double P, const double* be) -> EquilibriumResult
{
    return pimpl->solve(state, T, P, be);
}

auto EquilibriumSolver::solve(ChemicalState& state) -> EquilibriumResult
{
    return solve(state, state.temperature(), state.pressure(), state.elementAmounts());
}

auto EquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem, phaseIdentificationMethod phaseid, double _quantity) -> EquilibriumResult
{
	quantity = _quantity;
	phaseidMethod = phaseid;
    return solve(state, problem.temperature(), problem.pressure(), problem.elementAmounts());
}

auto EquilibriumSolver::properties() const -> const ChemicalProperties&
{
    return pimpl->properties;
}

auto EquilibriumSolver::sensitivity() -> const EquilibriumSensitivity&
{
    return pimpl->sensitivity();
}

auto EquilibriumSolver::dndT() -> VectorConstRef
{
    return pimpl->dndT();
}

auto EquilibriumSolver::dndP() -> VectorConstRef
{
    return pimpl->dndP();
}

auto EquilibriumSolver::dndb() -> VectorConstRef
{
    return pimpl->dndb();
}

} // namespace Reaktoro
