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

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

// Reactive transport test includes
#include <demos/cpp/DemosUtils.h>

using namespace Reaktoro;

/// Forward declaration
auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportResults& results) -> void;

int main()
{
    tic(TOTAL_TIME);

    // Step 1: Initialise auxiliary time-related constants
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;

    // Step 2: Define parameters for the reactive transport simulation
    ReactiveTransportParams params;

    // Define discretization parameters
    params.xl = 0.0;                                        // the x-coordinates of the left boundaries
    params.xr = 25.0;                                       // the x-coordinates of the right boundaries
    params.ncells = 243;                                    // the number of cells in the spacial discretization
    params.dx = (params.xr - params.xl) / params.ncells;    // the time step (in units of s)
    params.dt = 1*hour;                                     // the time step (in units of s)
    params.nsteps = 100;                                   // the total number of steps in the reactive transport simulation

    // Define physical and chemical parameters
    params.D = 0.0;             // the diffusion coefficient (in units of m2/s)
    params.v = 0.8e-5;          // the Darcy velocity (in units of m/s)
    params.T = 60.0;            // the temperature (in units of degC)
    params.P = 200 * 1.01325;   // the pressure (in units of bar)

    // Define the activity model for the aqueous species
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::HKF;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::Pitzer;
    params.activity_model = ReactiveTransportParams::AqueousActivityModel::DebyeHuckel;

    // Define activity model depending on the parameter
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // Define smart algorithm and related tolerances
    // -----------------------------------------------

    // Run smart algorithm with clustering
    params.method = SmartEquilibriumStrategy::Clustering;;
    params.smart_equilibrium_reltol = 1e-3;

//    // Run smart algorithm with priority queue
//    params.method = SmartEquilibriumStrategy::PriorityQueue;
//    params.smart_equilibrium_reltol = 2e-3;

//    // Run smart algorithm with nn search algorithm
//    params.method =  SmartEquilibriumStrategy::NearestNeighbour;
//    params.smart_equilibrium_reltol = 1e-1;

    // Output
    params.outputConsole();

    // Results
    ReactiveTransportResults results;

    // ----------------------------------------------------------- //
    // Execute reactive transport with different solvers
    // ----------------------------------------------------------- //

    // Smart equilibrium solver
    params.use_smart_equilibrium_solver = true;
    params.outputEquilibriumMethod();
    runReactiveTransport(params, results);

//    // Conventional equilibrium solver
//    params.outputEquilibriumMethod();
//    params.use_smart_equilibrium_solver = false;
//    runReactiveTransport(params, results);

    // Collect the time spent for total simulation (excluding search and store procedures costs)
    results.conventional_total = results.equilibrium_timing.solve;
    results.smart_total = results.smart_equilibrium_timing.solve;
    results.smart_total_ideal_search = results.smart_equilibrium_timing.solve
                                       - results.smart_equilibrium_timing.estimate_search
                                       - results.smart_equilibrium_timing.estimate_database_priority_update;
    results.smart_total_ideal_search_store = results.smart_equilibrium_timing.solve
                                             - results.smart_equilibrium_timing.estimate_search
                                             - results.smart_equilibrium_timing.estimate_database_priority_update
                                             - results.smart_equilibrium_timing.learn_storage;

    // Output speed-us
    std::cout << "speed up                            : "
              << results.conventional_total / results.smart_total << std::endl;
    std::cout << "speed up (with ideal search)        : "
              << results.conventional_total / results.smart_total_ideal_search << std::endl;
    std::cout << "speed up (with ideal search & store): "
              << results.conventional_total / results.smart_total_ideal_search_store << std::endl << std::endl;
    // Output reactive transport times and speedup
    std::cout << "time_reactive_transport_conventional: " << results.time_reactive_transport_conventional << std::endl;
    std::cout << "time_reactive_transport_smart       : " << results.time_reactive_transport_smart << std::endl;
    std::cout << "reactive_transport_speedup          : " << results.time_reactive_transport_conventional / results.time_reactive_transport_smart << std::endl;
    // Output total time
    std::cout << "total time                          : " << toc(TOTAL_TIME) << std::endl;

    return 0;
}

auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportResults& results) -> void
{
    // Step **: Create the results folder
    auto folder = params.makeResultsFolder("scaling-only-seawater");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;
    smart_equilibrium_options.amount_fraction_cutoff = params.amount_fraction_cutoff;
    smart_equilibrium_options.mole_fraction_cutoff = params.mole_fraction_cutoff;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor(database);

    // Define the list of selected elements
    StringList selected_elements = "H Cl S O Ba Ca Sr Na K Mg C Si";

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model == ReactiveTransportParams::AqueousActivityModel::HKF){
        // HKF full system
        editor.addAqueousPhaseWithElements(selected_elements);
    }
    else if(params.activity_model == ReactiveTransportParams::AqueousActivityModel::DebyeHuckel){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel(dhModel);
    }
    else if(params.activity_model == ActivityModel::Pitzer){
        // Pitzer full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    editor.addMineralPhase("Barite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // ************************************************************************************************************** //
    // Initial condition (IC)
    // ************************************************************************************************************** //

    // Define the initial condition *with formation water)
    EquilibriumInverseProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "atm");
    double water_kg = 1.0;
    problem_ic.add("H2O", water_kg, "kg");
    problem_ic.add("SO4", 10 * water_kg, "ug");
    problem_ic.add("Ca", 995 * water_kg, "mg");
    problem_ic.add("Ba", 995 * water_kg, "mg");
    problem_ic.add("Sr", 105 * water_kg, "mg");
    problem_ic.add("Na", 27250 * water_kg, "mg");
    problem_ic.add("K", 1730 * water_kg, "mg");
    problem_ic.add("Mg", 110 * water_kg, "mg");
    problem_ic.add("Cl", 45150 * water_kg, "mg");
    problem_ic.add("HCO3", 1980 * water_kg, "mg");
    problem_ic.pH(7.0, "HCl", "NaOH");

    // Equilibrate the initial condition
    ChemicalState state_ic = equilibrate(problem_ic);
    // Scale the percentage of the aqueous phase and the whole volume
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");

    // Define function to evaluate ph of the chemical system
    auto evaluate_pH = ChemicalProperty::pH(system);

    // Fetch the pH of the initial sate
    ChemicalProperties props = state_ic.properties();
    ChemicalScalar pH_ic = evaluate_pH(props);
    std::cout << "ph(FW) = " << pH_ic.val << std:: endl;

    // ************************************************************************************************************** //
    // Boundary condition (BC)
    // ************************************************************************************************************** //

    // Define the first boundary condition (with completion brine)
    EquilibriumProblem problem_bc_cb(system);
    problem_bc_cb.setTemperature(params.T, "celsius");
    problem_bc_cb.setPressure(params.P, "atm");
    problem_bc_cb.add("H2O", water_kg, "kg");
    problem_bc_cb.add("NaCl", 7, "mol");

    // Equilibrate the initial condition
    ChemicalState state_bc_cb = equilibrate(problem_bc_cb);
    // Scale the percentage of the aqueous phase and the whole volume
    state_bc_cb.scaleVolume(1.0, "m3");

    // Fetch the pH of the initial sate
    props = state_bc_cb.properties();
    ChemicalScalar pH_bc = evaluate_pH(props);
    std::cout << "ph(CB) = " << pH_bc.val << std:: endl;

    // Define the first boundary condition (with seawater)
    EquilibriumInverseProblem problem_bc_sw(system);
    problem_bc_sw.setTemperature(params.T, "celsius");
    problem_bc_sw.setPressure(params.P, "atm");
    problem_bc_sw.add("H2O", water_kg, "kg");
    problem_bc_sw.add("SO4--", 2710 * water_kg, "mg");
    problem_bc_sw.add("Ca++", 411 * water_kg, "mg");
    problem_bc_sw.add("Ba++", 0.01 * water_kg, "mg");
    problem_bc_sw.add("Sr++", 8 * water_kg, "mg");
    problem_bc_sw.add("Na+", 10760 * water_kg, "mg");
    problem_bc_sw.add("K+", 399 * water_kg, "mg");
    problem_bc_sw.add("Mg++", 1290 * water_kg, "mg");
    problem_bc_sw.add("Cl-", 19350 * water_kg, "mg");
    problem_bc_sw.add("HCO3-", 142 * water_kg, "mg");
    problem_bc_sw.pH(8.1, "HCl", "NaOH");

    // Equilibrate the initial condition
    ChemicalState state_bc_sw = equilibrate(problem_bc_sw);
    // Scale the percentage of the aqueous phase and the whole volume
    state_bc_sw.scaleVolume(1.0, "m3");

    // Fetch the pH of the initial sate
    props = state_bc_sw.properties();
    ChemicalScalar pH_sw = evaluate_pH(props);
    std::cout << "ph(SW) = " << pH_sw.val << std:: endl;

    // Step **: Create the mesh for the column
    Mesh mesh(params.ncells, params.xl, params.xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step **: Define the options for the reactive transport solver
    ReactiveTransportOptions reactive_transport_options;
    reactive_transport_options.use_smart_equilibrium_solver = params.use_smart_equilibrium_solver;
    reactive_transport_options.equilibrium = equilibrium_options;
    reactive_transport_options.smart_equilibrium = smart_equilibrium_options;

    // Step **: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(system);
    rtsolver.setOptions(reactive_transport_options);
    rtsolver.setMesh(mesh);
    rtsolver.setVelocity(params.v);
    rtsolver.setDiffusionCoeff(params.D);
    rtsolver.setBoundaryState(state_bc_sw);
    rtsolver.setTimeStep(params.dt);
    rtsolver.initialize();

    // Step **: Define the quantities that should be output for every cell, every time step
    ChemicalOutput output(rtsolver.output());
    output.add("pH");
    output.add("speciesMolality(H+)");
    output.add("speciesMolality(Cl-)");
    output.add("speciesMolality(SO4--)");
    output.add("speciesMolality(Ba++)");
    output.add("speciesMolality(Ca++)");
    output.add("speciesMolality(Sr++)");
    output.add("speciesMolality(Na+)");
    output.add("speciesMolality(Barite)");
    output.add("elementmolality(Ba)");
    output.add("elementmolality(C)");
    output.add("elementmolality(Ca)");
    output.add("elementmolality(Cl)");
    output.add("elementmolality(H)");
    output.add("elementmolality(K)");
    output.add("elementmolality(Mg)");
    output.add("elementmolality(Na)");
    output.add("elementmolality(O)");
    output.add("elementmolality(S)");
    output.add("elementmolality(Si)");
    output.add("elementmolality(Sr)");
    output.add("elementmolality(Z)");
    output.add("phaseAmount(Barite)");
    output.add("phaseMass(Barite)");
    output.add("phaseVolume(Barite)");
    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(REACTIVE_TRANSPORT_STEPS);

    // Reactive transport simulations in the cycle
    while (step <  params.nsteps)
    {
        // Print the progress of simulations
        //std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }

    // Print the content of the cluster if the smart equilibrium with clustering is used
    if(params.use_smart_equilibrium_solver)
        rtsolver.outputSmartSolverInfo();

    // Stop the time for the reactive transport simulation
    if(params.use_smart_equilibrium_solver) results.time_reactive_transport_smart = toc(REACTIVE_TRANSPORT_STEPS);
    else results.time_reactive_transport_conventional = toc(REACTIVE_TRANSPORT_STEPS);

    // Step **: Collect the analytics related to reactive transport performance
    auto analysis = profiler.analysis();
    auto rt_results = profiler.results();

    // Step **: Generate json output file with collected profiling data
    if(params.use_smart_equilibrium_solver)  JsonOutput(folder + "/" + "analysis-smart.json") << analysis;
    else    JsonOutput(folder + "/" + "analysis-conventional.json") << analysis;

    // Step **: Save equilibrium timing to compare the speedup of smart equilibrium solver versus conventional one
    if(params.use_smart_equilibrium_solver) {
        results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
        results.smart_equilibrium_acceptance_rate = analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;

        std::cout << "smart equilibrium acceptance rate   : " << results.smart_equilibrium_acceptance_rate << " / "
                  << (1 - results.smart_equilibrium_acceptance_rate) * params.ncells *params.nsteps
                  << " fully evaluated GEMS out of " << params.ncells * params.nsteps  << std::endl;

    }
    else results.equilibrium_timing = analysis.equilibrium.timing;
}
