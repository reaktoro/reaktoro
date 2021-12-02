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

// C++ includes
#include<fstream>

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
    int week = 7 * day;

    // Step 2: Define parameters for the reactive transport simulation
    ReactiveTransportParams params;

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 1.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    params.nsteps = 1; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 30 * minute; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 1.0e-9;     // the diffusion coefficient (in units of m2/s)
    params.v = 1.0 / week; // the Darcy velocity (in units of m/s)
    params.T = 60.0;       // the temperature (in units of degC)
    params.P = 100;        // the pressure (in units of bar)

    // Define the activity model for the aqueous species
    params.activity_model = RTActivityModel::HKF;
    //params.activity_model = RTActivityModel::Pitzer;
    //params.activity_model = RTActivityModel::DebyeHuckel;

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


    std::ifstream readfile;
    readfile.open("rt-dolomitization-dt-1800-ncells-100-nsteps-1-hkf-conventional/step-0.txt");
    std::string output_string;
    char output[100];

    if(readfile.is_open())
    {
//        readfile >> output;
//        std::cout << output << std::endl;
//        readfile >> output;
//        std::cout << output << std::endl;

        getline(readfile, output_string);
        std::cout << output_string << std::endl;

        readfile >> output;
        std::cout << output << std::endl;
        readfile >> output;
        std::cout << output << std::endl;
        readfile >> output;
        std::cout << output << std::endl;
        readfile >> output;
        std::cout << output << std::endl;
        readfile >> output;
        std::cout << output << std::endl;
        readfile >> output;
        std::cout << output << std::endl;
        readfile >> output;
        std::cout << output << std::endl;
        readfile >> output;
        std::cout << output << std::endl;
        readfile >> output;
        std::cout << output << std::endl;
    }

    return 0;
}

auto runReactiveTransport(ReactiveTransportParams& params, ReactiveTransportResults& results) -> void
{
    // Step **: Create the results folder
    auto folder = params.makeResultsFolder("dolomitization");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;
    equilibrium_options.hessian = params.hessian;

    // Define the list of selected species
    StringList selected_species = "H2O(aq) H+ OH- Na+ Cl- NaCl(aq) Ca+2 Mg+2 HCO3- CO2(aq) CO3-2 CaCl+ CaCl2(aq) O2(aq) Ca(HCO3)+ MgCl+ Mg(HCO3)+";

    // Initialize a thermodynamic database
    SupcrtDatabase db("supcrtbl");

    // Create an aqueous phase automatically selecting all species with given elements, excluding species with tag `organic`
    AqueousPhase aqueousphase(selected_species);

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model == RTActivityModel::HKF)
    {
        // HKF full system
        aqueousphase.setActivityModel(chain(
                ActivityModelHKF(),
                ActivityModelDrummond("CO2")
        ));
    }
    else if(params.activity_model == RTActivityModel::Pitzer)
    {
        // Pitzer full system
        aqueousphase.setActivityModel(chain(
                ActivityModelPitzerHMW(),
                ActivityModelDrummond("CO2")
        ));
    }
    else if(params.activity_model == RTActivityModel::DebyeHuckel)
    {
        // Pitzer full system
        aqueousphase.setActivityModel(chain(
                ActivityModelDebyeHuckel(),
                ActivityModelDrummond("CO2")
        ));
    }

    // Create a mineral phase
    MineralPhases mineralphases("Quartz Calcite Dolomite");

    // Collect all above-defined phases
    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(mineralphases);

    // Construct the chemical system
    ChemicalSystem system(phases);

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);

    // Step **: Create the Partition of inert and equilibrium species
    Partition partition(system);

    // Define aqueous and chemical properties
    AqueousProps aprops(system);
    ChemicalProps props(system);

    // Define initial equilibrium state
    ChemicalState state_ic(system);
    state_ic.setTemperature(params.T, "celsius");
    state_ic.setPressure(params.P, "bar");
    state_ic.setSpeciesMass("H2O(aq)", 1.00, "kg");
    state_ic.setSpeciesAmount("O2(aq)", 1.00, "umol");
    state_ic.setSpeciesAmount("NaCl(aq)", 0.70, "mol");
    state_ic.setSpeciesAmount("MgCl+", 1e-10, "mol");
    state_ic.setSpeciesAmount("Cl-"  , 1e-10, "mol");
    state_ic.setSpeciesAmount("Calcite", 10.00, "mol");
    state_ic.setSpeciesAmount("Quartz", 10.00, "mol");

    solver.solve(state_ic);
    aprops.update(state_ic);
    std::cout << "pH (IC) = " << aprops.pH() << std::endl;

    // Define boundary equilibrium state
    ChemicalState state_bc(system);
    state_bc.setTemperature(params.T, "celsius");
    state_bc.setPressure(params.P, "bar");
    state_bc.setSpeciesMass("H2O(aq)", 1.00, "kg");
    state_bc.setSpeciesAmount("O2(aq)", 1.00, "umol");
    state_bc.setSpeciesAmount("NaCl(aq)", 0.90, "mol");
    state_bc.setSpeciesAmount("MgCl+", 0.05, "mol");
    state_bc.setSpeciesAmount("Cl-"  , 0.05, "mol");
    state_bc.setSpeciesAmount("CaCl2(aq)", 0.01, "mol");
    state_bc.setSpeciesAmount("CO2(aq)", 0.75, "mol");

    solver.solve(state_bc);
    aprops.update(state_bc);
    std::cout << "pH (BC) = " << aprops.pH() << std::endl;

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("AqueousPhase", 0.1, "m3");    // 10% if the 1.0m3
    state_ic.scalePhaseVolume("Quartz", 0.882, "m3");   // 0.882 = 0.98 * 0.9 (0.9 is 90% of 1.0m3, 0.98 is 98% quartz of the rock)
    state_ic.scalePhaseVolume("Calcite", 0.018, "m3");  // 0.018 = 0.02 * 0.9 (0.9 is 90% of 1.0m3, 0.02 is 2% calcite of the rock)

    // Step **: Create the mesh for the column
    Mesh mesh(params.ncells, params.xl, params.xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step **: Define the options for the reactive transport solver
    ReactiveTransportOptions reactive_transport_options;
    reactive_transport_options.use_smart_solver = params.use_smart_equilibrium_solver;
    reactive_transport_options.equilibrium = equilibrium_options;

    // Step **: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(partition);
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
    output.add("speciesMolality(H+)");
    output.add("speciesMolality(Ca+2)");
    output.add("speciesMolality(Mg+2)");
    output.add("speciesMolality(HCO3-)");
    output.add("speciesMolality(CO2(aq))");
    output.add("speciesMolality(CO3-2)");
    output.add("speciesMolality(CaCl+)");
    output.add("speciesMolality(Ca(HCO3)+)");
    output.add("speciesMolality(MgCl+)");
    output.add("speciesMolality(Mg(HCO3)+)");
    output.add("speciesMolality(OH-)");
    output.add("phaseVolume(Calcite)");
    output.add("phaseVolume(Dolomite)");
    output.add("elementMolality(C)");
    output.add("elementMolality(Ca)");
    output.add("elementMolality(Cl)");
    output.add("elementMolality(H)");
    output.add("elementMolality(Mg)");
    output.add("elementMolality(Na)");
    output.add("elementMolality(O)");
    output.add("elementMolality(Si)");
    output.add("speciesAmount(H+)");
    output.add("speciesAmount(Ca+2)");
    output.add("speciesAmount(Mg+2)");
    output.add("speciesAmount(HCO3-)");
    output.add("speciesAmount(CO2(aq))");
    output.add("speciesAmount(CO3-2)");
    output.add("speciesAmount(CaCl+)");
    output.add("speciesAmount(Ca(HCO3)+)");
    output.add("speciesAmount(MgCl+)");
    output.add("speciesAmount(Mg(HCO3)+)");
    output.add("speciesAmount(OH-)");
    output.add("speciesAmount(Calcite)");
    output.add("speciesAmount(Dolomite)");
    output.filename(folder + "/" + "step.txt");

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
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }
    std::cout << std::endl;

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