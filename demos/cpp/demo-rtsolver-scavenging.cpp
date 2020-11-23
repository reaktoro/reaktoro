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
#include <demos/cpp/TestUtils.h>

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
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 100.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    params.nsteps = 100; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 0.05*day; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 0.0;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.05e-5; // the Darcy velocity (in units of m/s)
    params.T = 25.0;                     // the temperature (in units of degC)
    params.P = 1.01325;                      // the pressure (in units of bar)

    // Define parameters of the equilibrium solvers
    params.smart_equilibrium_reltol = 0.01;

    // Define the activity model for the aqueous species
    params.activity_model = ReactiveTransportParams::AqueousActivityModel::HKF;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::HKFSelectedSpecies;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::Pitzer;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::PitzerSelectedSpecies;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::DebyeHuckel;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::DebyeHuckelSelectedSpecies;

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

    // Define equilibrium solver cutoff tolerances
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // Output
    params.outputConsole();

    // Results
    ReactiveTransportResults results;

    // Execute reactive transport with different solvers
    params.use_smart_equilibrium_solver = true; runReactiveTransport(params, results);
    params.use_smart_equilibrium_solver = false; runReactiveTransport(params, results);

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
    auto folder = params.makeResultsFolder("scavenging");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;
    smart_equilibrium_options.amount_fraction_cutoff = params.amount_fraction_cutoff;
    smart_equilibrium_options.mole_fraction_cutoff = params.mole_fraction_cutoff;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    Database database("supcrt07.xml");

    // Step **: Define Debye-Huckel parameters
    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor(database);

    // Define the list of selected species
    StringList selected_species = {"Ca(HCO3)+", "CO3--", "CaCO3(aq)", "Ca++", "CaSO4(aq)", "CaOH+", "Cl-",
                                   "FeCl++", "FeCl2(aq)", "FeCl+", "Fe++", "FeOH+",  "FeOH++", "Fe+++",
                                   "H2(aq)", "HSO4-", "H2S(aq)", "HS-", "H2O(l)",  "H+", "OH-", "HCO3-",
                                   "K+", "KSO4-",
                                   "Mg++", "MgSO4(aq)", "MgCO3(aq)", "MgOH+", "Mg(HCO3)+",
                                   "Na+", "NaSO4-",
                                   "O2(aq)",
                                   "S5--", "S4--", "S3--", "S2--", "SO4--"};
    // Define the list of selected elements
    StringList selected_elements = "C Ca Cl Fe H K Mg Na O S";

    // Define activity model depending on the parameter
    if(params.activity_model == ReactiveTransportParams::AqueousActivityModel::HKF){
        // HKF full system
        editor.addAqueousPhaseWithElements(selected_elements);
    }
    else if(params.activity_model == ReactiveTransportParams::AqueousActivityModel::HKFSelectedSpecies){
        // HKF with selected species
        editor.addAqueousPhase(selected_species);
    }
    else if(params.activity_model == ReactiveTransportParams::AqueousActivityModel::PitzerSelectedSpecies){
        // Pitzer full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == ReactiveTransportParams::AqueousActivityModel::PitzerSelectedSpecies){
        // Pitzer with selected species
        editor.addAqueousPhase(selected_species)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == ReactiveTransportParams::AqueousActivityModel::DebyeHuckelSelectedSpecies){
        // Debye-Huckel with selected species
        editor.addAqueousPhase(selected_species)
                .setChemicalModelDebyeHuckel(dhModel);
    }
    else if(params.activity_model == ReactiveTransportParams::AqueousActivityModel::DebyeHuckel){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel(dhModel);
    }
    editor.addMineralPhase("Pyrrhotite"); // also known as Mackinawite
    editor.addMineralPhase("Siderite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumInverseProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O", 58.0, "kg");
    problem_ic.add("Cl-", 1122.3e-3, "kg");
    problem_ic.add("Na+", 624.08e-3, "kg");
    problem_ic.add("SO4--", 157.18e-3, "kg");
    problem_ic.add("Mg++", 74.820e-3, "kg");
    problem_ic.add("Ca++", 23.838e-3, "kg");
    problem_ic.add("K+", 23.142e-3, "kg");
    problem_ic.add("HCO3-", 8.236e-3, "kg");
    problem_ic.add("O2(aq)", 58e-12, "kg");
    problem_ic.add("Pyrrhotite", 0.0, "mol");
    problem_ic.add("Siderite", 0.5, "mol");
    problem_ic.pH(8.951);
    problem_ic.pE(8.676);

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumInverseProblem problem_bc(system);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O", 58.0, "kg");
    problem_bc.add("Cl-", 1122.3e-3, "kg");
    problem_bc.add("Na+", 624.08e-3, "kg");
    problem_bc.add("SO4--", 157.18e-3, "kg");
    problem_bc.add("Mg++", 74.820e-3, "kg");
    problem_bc.add("Ca++", 23.838e-3, "kg");
    problem_bc.add("K+", 23.142e-3, "kg");
    problem_bc.add("HCO3-", 8.236e-3, "kg");
    problem_bc.add("O2(aq)", 58e-12, "kg");
    problem_bc.add("Pyrrhotite", 0.0, "mol");
    problem_bc.add("Siderite", 0.0, "mol");
    problem_bc.add("HS-", 0.0196504, "mol");
    problem_bc.add("H2S(aq)", 0.167794, "mol");
    problem_bc.pH(5.726);
    problem_bc.pE(8.220);

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");

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
    rtsolver.setBoundaryState(state_bc);
    rtsolver.setTimeStep(params.dt);
    rtsolver.initialize();

    // Step **: Define the quantities that should be output for every cell, every time step
    ChemicalOutput output(rtsolver.output());
    output.add("pH");
    output.add("speciesAmount(H+)");
    output.add("speciesAmount(HS-)");
    output.add("speciesAmount(S2--)");
    output.add("speciesAmount(CO3--)");
    output.add("speciesAmount(HSO4-)");
    output.add("speciesAmount(H2S(aq))");
    output.add("phaseAmount(Pyrrhotite)");
    output.add("phaseAmount(Siderite)");
    output.add("phaseVolume(Pyrrhotite)");
    output.add("phaseVolume(Siderite)");
    output.add("elementAmount(C)");
    output.add("elementAmount(Ca)");
    output.add("elementAmount(Cl)");
    output.add("elementAmount(Fe)");
    output.add("elementAmount(H)");
    output.add("elementAmount(K)");
    output.add("elementAmount(Mg)");
    output.add("elementAmount(Na)");
    output.add("elementAmount(O)");
    output.add("elementAmount(S)");
    output.add("elementAmount(Z)");
    output.add("speciesAmount(Fe++)");
    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(REACTIVE_TRANSPORT_STEPS);

    // Reactive transport simulations in the cycle
    while (step < params.nsteps)
    {
        // Print the progress of simulations
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

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

