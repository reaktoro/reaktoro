// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Svetlana Kyas (2 December, 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

// Reactive transport test includes
#include <examples/resources/DemoUtils.hpp>

using namespace Reaktoro;

/// Forward declaration
auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportResults& results) -> void;
auto outputStatistics(ReactiveTransportResults& results) -> void;

int main()
{
    tic(TOTAL_TIME)

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
    params.nsteps = 1000; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 0.05*day; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 0.0;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.05e-5; // the Darcy velocity (in units of m/s)
    params.T = 25.0;                     // the temperature (in units of degC)
    params.P = 1.01325;                      // the pressure (in units of bar)

    // Define the activity model for the aqueous species
    //params.activity_model = RTActivityModel::HKF;
    //params.activity_model = RTActivityModel::Pitzer;
    params.activity_model = RTActivityModel::DebyeHuckel;

    // ----------------------------------------------------------- //
    // Optimization algorithm parameters
    // ----------------------------------------------------------- //

    params.hessian = GibbsHessian::Exact;
    //params.hessian = GibbsHessian::Approximation;

    // Output
    params.outputConsole();

    // Results
    ReactiveTransportResults results;

    // ----------------------------------------------------------- //
    // Execute reactive transport with different solvers
    // ----------------------------------------------------------- //

    // Conventional equilibrium solver
    params.outputEquilibriumMethod();
    params.use_smart_equilibrium_solver = false;
    runReactiveTransport(params, results);

    outputStatistics(results);

    // Output total time
    std::cout << "total time                          : " << toc(TOTAL_TIME) << std::endl;

    return 0;
}

auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportResults& results) -> void
{
    // Step **: Create the results folder
    auto folder = params.makeResultsFolder("scavenging-with-hematite");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;
    equilibrium_options.hessian = params.hessian;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    SupcrtDatabase db("supcrt07");

    // Define the list of selected species
    StringList selected_species = {"Ca(HCO3)+", "CO3-2", "CaCO3(aq)", "Ca+2", "CaSO4(aq)", "CaOH+", "Cl-",
                                   "FeCl+2", "FeCl2(aq)", "FeCl+", "Fe+2", "FeOH+",  "FeOH+2", "Fe+3",
                                   "H2(aq)", "HSO4-", "H2S(aq)", "HS-", "H2O(aq)",  "H+", "OH-", "HCO3-",
                                   "K+", "KSO4-", "Mg+2", "MgSO4(aq)", "MgCO3(aq)", "MgOH+", "Mg(HCO3)+",
                                   "Na+", "NaSO4-", "O2(aq)", "S5-2", "S4-2", "S3-2", "S2-2", "SO4-2"};

    // Create an aqueous phase automatically selecting all species with given elements, excluding species with tag `organic`
    AqueousPhase aqueousphase(selected_species);

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model == RTActivityModel::HKF)
    {
        // HKF full system
        aqueousphase.setActivityModel(ActivityModelHKF());
    }
    else if(params.activity_model == RTActivityModel::DebyeHuckel)
    {
        // Debye-Huckel full system
        aqueousphase.setActivityModel(ActivityModelDebyeHuckelPHREEQC());
    }
    else if(params.activity_model == RTActivityModel::Pitzer)
    {
        // Pitzer full system
        aqueousphase.setActivityModel(ActivityModelPitzerHMW());
    }

    // Create a mineral phase
    MineralPhases mineralphases("Pyrite Hematite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(db, aqueousphase, mineralphases);

    // Define equilibrium specs
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    // Define aqueous and chemical properties
    AqueousProps aprops(system);
    ChemicalProps props(system);

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(specs);

    // ************************************************************************************************************** //
    // Initial condition (IC)
    // ************************************************************************************************************** //

    EquilibriumConditions conditions_ic(specs);
    conditions_ic.temperature(params.T, "celsius");
    conditions_ic.pressure(params.P, "bar");
    conditions_ic.pH(8.951); // H+ will be added as much as needed to achieve this pH

    // Define the initial condition (with formation water)
    ChemicalState state_ic(system);
    state_ic.setTemperature(params.T, "celsius");
    state_ic.setPressure(params.P, "bar");
    state_ic.setSpeciesMass("H2O(aq)", 58.0, "kg");
    state_ic.setSpeciesMass("Cl-", 1122.3e-3, "kg");
    state_ic.setSpeciesMass("Na+", 624.08e-3, "kg");
    state_ic.setSpeciesMass("SO4-2", 157.18e-3, "kg");
    state_ic.setSpeciesMass("Mg+2", 74.820e-3, "kg");
    state_ic.setSpeciesMass("Ca+2", 23.838e-3, "kg");
    state_ic.setSpeciesMass("K+", 23.142e-3, "kg");
    state_ic.setSpeciesMass("HCO3-", 8.236e-3, "kg");
    state_ic.setSpeciesMass("O2(aq)", 58e-12, "kg");
    state_ic.setSpeciesAmount("Pyrite", 0.0, "mol");
    state_ic.setSpeciesAmount("Hematite", 0.5, "mol");

    solver.solve(state_ic, conditions_ic);
    aprops.update(state_ic);
    std::cout << "pH (IC) = " << aprops.pH() << std::endl;
    state_ic.scalePhaseVolume("AqueousPhase", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scaleVolume(1.0, "m3");

    // ************************************************************************************************************** //
    // Boundary condition (BC)
    // ************************************************************************************************************** //

    EquilibriumConditions conditions_bc(specs);
    conditions_bc.temperature(params.T, "celsius");
    conditions_bc.pressure(params.P, "bar");
    conditions_bc.pH(5.726); // H+ will be added as much as needed to achieve this pH

    // Define the initial condition (with seawater water)
    ChemicalState state_bc(system);
    state_bc.setTemperature(params.T, "celsius");
    state_bc.setPressure(params.P, "bar");
    state_bc.setSpeciesMass("H2O(aq)", 58.0, "kg");
    state_bc.setSpeciesMass("Cl-", 1122.3e-3, "kg");
    state_bc.setSpeciesMass("Na+", 624.08e-3, "kg");
    state_bc.setSpeciesMass("SO4-2", 157.18e-3, "kg");
    state_bc.setSpeciesMass("Mg+2", 74.820e-3, "kg");
    state_bc.setSpeciesMass("Ca+2", 23.838e-3, "kg");
    state_bc.setSpeciesMass("K+", 23.142e-3, "kg");
    state_bc.setSpeciesMass("HCO3-", 8.236e-3, "kg");
    state_bc.setSpeciesMass("O2(aq)", 58e-12, "kg");
    state_bc.setSpeciesAmount("Pyrite", 0.0, "mol");
    state_bc.setSpeciesAmount("Hematite", 0.0, "mol");
    state_bc.setSpeciesAmount("HS-", 0.0196504, "mol");
    state_bc.setSpeciesAmount("H2S(aq)", 0.167794, "mol");

    solver.solve(state_bc, conditions_bc);
    aprops.update(state_bc);
    std::cout << "pH (BC) = " << aprops.pH() << std::endl;
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Create the mesh for the column
    Mesh mesh(params.ncells, params.xl, params.xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step **: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(system);
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
    output.add("speciesAmount(S2-2)");
    output.add("speciesAmount(CO3-2)");
    output.add("speciesAmount(HSO4-)");
    output.add("speciesAmount(H2S(aq))");
    output.add("speciesAmount(Fe+2)");
    output.add("speciesAmount(Fe+3)");
    output.add("speciesAmount(Pyrite)");
    output.add("speciesAmount(Hematite)");
    output.add("phaseVolume(Pyrite)");
    output.add("phaseVolume(Hematite)");
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

    // Stop the time for the reactive transport simulation
    if(params.use_smart_equilibrium_solver) results.time_reactive_transport_smart = toc(REACTIVE_TRANSPORT_STEPS);
    else results.time_reactive_transport_conventional = toc(REACTIVE_TRANSPORT_STEPS);

    // Step **: Collect the analytics related to reactive transport performance
    auto analysis = profiler.analysis();
    auto rt_results = profiler.results();

    // Step **: Generate json output file with collected profiling data
    if(params.use_smart_equilibrium_solver)  JsonOutput(folder + "/" + "analysis-smart.json") << analysis;
    else    JsonOutput(folder + "/" + "analysis-conventional.json") << analysis;

    results.equilibrium_timing = analysis.equilibrium.timing;
}

auto outputStatistics(ReactiveTransportResults& results) -> void
{
    // Collect the time spent for total simulation (excluding search and store procedures costs)
    results.conventional_total = results.equilibrium_timing.solve;

    // Output reactive transport times and speedup
    std::cout << "time_reactive_transport_conventional: " << results.time_reactive_transport_conventional << std::endl;
}