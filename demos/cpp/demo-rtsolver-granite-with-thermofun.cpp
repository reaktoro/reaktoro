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
#include <ThermoFun/ThermoFun.h>

using namespace Reaktoro;

/// Forward declaration
auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportResults& results) -> void;

int main()
{
    tic(TOTAL_TIME)

    // Step 1: Initialise auxiliary time-related constants
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;
    int week = 7 * day;

    // Step 2: Define parameters for the reactive transport simulation
    ReactiveTransportParams params;

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 1.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    params.nsteps = 10000; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 30 * minute; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 1.0e-9;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.0 / week; // the Darcy velocity (in units of m/s)
    params.T = 300;        // the temperature (in units of degC)
    // Calculate the water saturation pressure using the Wagner and Pruss (1995) equation of state
    params.P = Reaktoro::waterSaturatedPressureWagnerPruss(Temperature(params.T + 273.15)).val * 1e-5;

    // Define the activity model for the aqueous species
    //params.activity_model = ActivityModel::HKF;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::HKF;
    params.activity_model = ReactiveTransportParams::AqueousActivityModel::HKFSelectedSpecies;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::Pitzer;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::PitzerSelectedSpecies;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::DebyeHuckel;
    //params.activity_model = ReactiveTransportParams::AqueousActivityModel::DebyeHuckelSelectedSpecies;

    // Define activity model depending on the parameter
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // Define smart algorithm and related tolerances
    // -----------------------------------------------

    // Run smart algorithm with clustering
    params.method = SmartEquilibriumStrategy::Clustering;
    params.smart_equilibrium_reltol = 1e-2;

//    // Run smart algorithm with priority queue
//    params.smart_method = SmartEquilibriumStrategy::PriorityQueue;
//    params.smart_equilibrium_reltol = 2e-3;

//    // Run smart algorithm with nn search algorithm
//    params.smart_method =  SmartEquilibriumStrategy::NearestNeighbour;
//    params.smart_equilibrium_reltol = 1e-1;

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
    std::cout << "smart equilibrium acceptance rate   : " << results.smart_equilibrium_acceptance_rate << " / "
              << (1 - results.smart_equilibrium_acceptance_rate) * params.ncells *params.nsteps
              << " fully evaluated GEMS out of " << params.ncells * params.nsteps  << std::endl;
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
    auto folder = params.makeResultsFolder("granite-thermofun");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;

    // Step **: Define smart chemical equilibrium solver options
    SmartEquilibriumOptions smart_equilibrium_options;
    smart_equilibrium_options.reltol = params.smart_equilibrium_reltol;
    smart_equilibrium_options.amount_fraction_cutoff = params.amount_fraction_cutoff;
    smart_equilibrium_options.mole_fraction_cutoff = params.mole_fraction_cutoff;
    smart_equilibrium_options.method = params.method;


    ThermoFun::Database database("databases/thermofun/aq17-thermofun.json");

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor(database);
    editor.setTemperatures({params.T}, "celsius");
    editor.setPressures({params.P}, "bar");

    // Define the list of selected elements
    StringList selected_elements = "Al Cl H K Na O Si";

    // Define the list of selected species
    StringList selected_species = "H2O@ H+ OH- Cl- HCl@ "
                                  "Na+ NaOH@ NaHSiO3@ NaCl@ NaAl(OH)4@ "
                                  "K+ KOH@ KCl@ KAlO2@ "
                                  "Al+3 AlOH+2 Al(OH)2+ Al(OH)3@ Al(OH)4-";

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model ==  ReactiveTransportParams::AqueousActivityModel::HKF){
        // HKF full system
        editor.addAqueousPhaseWithElements(selected_elements);
    }
    else if(params.activity_model ==  ReactiveTransportParams::AqueousActivityModel::HKFSelectedSpecies){
        // HKF selected species
        editor.addAqueousPhase(selected_species);
    }
    else if(params.activity_model ==  ReactiveTransportParams::AqueousActivityModel::Pitzer){
        // Pitzer full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model ==  ReactiveTransportParams::AqueousActivityModel::PitzerSelectedSpecies){
        // Pitzer selected species
        editor.addAqueousPhase(selected_species)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model ==  ReactiveTransportParams::AqueousActivityModel::DebyeHuckel){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model ==  ReactiveTransportParams::AqueousActivityModel::DebyeHuckelSelectedSpecies){
        // Debye-Huckel selected species
        editor.addAqueousPhase(selected_species)
                .setChemicalModelDebyeHuckel()
                .setActivityModelDrummondCO2();
    }
    editor.addMineralPhase("Quartz"); // SiO2
    editor.addMineralPhase("Diaspore"); // AlO(OH)
    editor.addMineralPhase("Gibbsite"); // Al(OH)3
    editor.addMineralPhase("Andalusite"); // Al2SiO5
    editor.addMineralPhase("Kyanite"); // Al2SiO5
    editor.addMineralPhase("Sillimanite"); // Al2SiO5
    editor.addMineralPhase("Muscovite"); // KAl2(AlSi3)O10(OH)2
    editor.addMineralPhase("Paragonite"); // NaAl2(AlSi3)O10(OH)2
    editor.addMineralPhase("Pyrophyllite"); // Al2Si4O10(OH)2
    editor.addMineralPhase("Kaolinite"); // Al2Si2O5(OH)4
    editor.addMineralPhase("Albite"); // Na(AlSi3)O8
    editor.addMineralPhase("Microcline"); // K(AlSi3)O8

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O", 55.51, "mol"); // H2O 55.51 M
    problem_ic.add("NaCl", 0.27, "mol"); // NaCl (aq) 0.27 M
    problem_ic.add("KCl", 0.03, "mol"); // KCl (aq)  0.03 M

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumProblem problem_bc(system);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O", 55.51, "mol"); // H2O 55.51 M
    problem_bc.add("HCl", 0.1, "mol"); // HCl (aq) 0.1 M
    problem_bc.add("NaCl", 0.17, "mol"); // NaCl (aq) 0.17 M
    problem_bc.add("KCl", 0.03, "mol"); // KCl (aq)  0.03 M

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    // Define function to evaluate ph of the chemical system
    auto evaluate_pH = ChemicalProperty::pH(system);

    // Fetch the pH of the initial sate
    ChemicalProperties props = state_ic.properties();
    ChemicalScalar pH_ic = evaluate_pH(props);
    std::cout << "ph(IC) = " << pH_ic.val << std:: endl;

    // Fetch the pH of the boundary sate
    ChemicalProperties props_bc = state_bc.properties();
    ChemicalScalar pH_bc = evaluate_pH(props_bc);
    std::cout << "ph(BC) = " << pH_bc.val << std:: endl;

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scalePhaseVolume("Quartz", 0.3 * 0.9, "m3");    // 30% of 90% of remaining volume
    state_ic.scalePhaseVolume("Muscovite", 0.05 * 0.9, "m3");    // 5% of 90% of remaining volume
    state_ic.scalePhaseVolume("Albite", 0.33 * 0.9, "m3");    // 33% of 90% of remaining volume
    state_ic.scalePhaseVolume("Microcline", 0.32 * 0.9, "m3");    // 32% of 90% of remaining volume

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
    output.add("speciesMolality(Cl-)");
    output.add("speciesMolality(HCl@)");
    output.add("speciesMolality(Na+)");
    output.add("speciesMolality(NaCl@)");
    output.add("speciesMolality(OH-)");
    output.add("speciesMolality(NaOH@)");
    output.add("speciesMolality(NaHSiO3@)");
    output.add("speciesMolality(K+)");
    output.add("speciesMolality(KOH@)");
    output.add("speciesMolality(KCl@)");
    output.add("speciesMolality(Al+3)");
    output.add("speciesMolality(AlOH+2)");
    output.add("speciesMolality(KAlO2@)");
    output.add("speciesMolality(NaAl(OH)4@)");
    output.add("speciesMolality(Al(OH)2+)");
    output.add("speciesMolality(Al(OH)3@)");
    output.add("speciesMolality(Al(OH)4-)");
    output.add("speciesMolality(Quartz)");
    output.add("speciesMolality(Diaspore)");
    output.add("speciesMolality(Gibbsite)");
    output.add("speciesMolality(Andalusite)");
    output.add("speciesMolality(Kyanite)");
    output.add("speciesMolality(Sillimanite)");
    output.add("speciesMolality(Muscovite)");
    output.add("speciesMolality(Paragonite)");
    output.add("speciesMolality(Pyrophyllite)");
    output.add("speciesMolality(Kaolinite)");
    output.add("speciesMolality(Albite)");
    output.add("speciesMolality(Microcline)");
    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(REACTIVE_TRANSPORT_STEPS)

     // Reactive transport simulations in the cycle
    while (step < params.nsteps)
    {
        // Print the progress of calculations
        //std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }
    std::cout << std::endl;

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

