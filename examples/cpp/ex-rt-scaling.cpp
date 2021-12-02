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

    // Step 2: Define parameters for the reactive transport simulation
    ReactiveTransportParams params;

    // Define discretization parameters
    params.xl = 0.0;                                        // the x-coordinates of the left boundaries
    params.xr = 25.0;                                       // the x-coordinates of the right boundaries
    params.ncells = 243;                                    // the number of cells in the spacial discretization
    params.dx = (params.xr - params.xl) / params.ncells;    // the time step (in units of s)
    params.dt = 1*hour;                                     // the time step (in units of s)
    params.nsteps = 200;                                   // the total number of steps in the reactive transport simulation

    // Define physical and chemical parameters
    params.D = 0.0;             // the diffusion coefficient (in units of m2/s)
    params.v = 0.8e-5;          // the Darcy velocity (in units of m/s)
    params.T = 60.0;            // the temperature (in units of degC)
    params.P = 200 * 1.01325;   // the pressure (in units of bar)

    // Define the activity model for the aqueous species
    //params.activity_model = RTActivityModel::HKF;
    //params.activity_model = RTActivityModel::Pitzer;
    params.activity_model = RTActivityModel::DebyeHuckel;

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
    auto folder = params.makeResultsFolder("scaling-only-seawater");

    // Step **: Define chemical equilibrium solver options
    EquilibriumOptions equilibrium_options;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    SupcrtDatabase db("supcrtbl");

    // Define the list of selected species
    StringList selected_species = "H2O(aq) H+ OH- Na+ Cl- NaCl(aq) Ca+2 Mg+2 Ba+2 Sr+2 K+ HCO3- CO2(aq) CO3-2 CaCl+ CaCl2(aq) O2(aq) Ca(HCO3)+ MgCl+ Mg(HCO3)+ SO4-2";

    // Create an aqueous phase automatically selecting all species with given elements, excluding species with tag `organic`
    AqueousPhase aqueousphase(selected_species);

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model == RTActivityModel::HKF)
    {
        // HKF full system
        aqueousphase.setActivityModel(chain(
                ActivityModelHKF(),
                ActivityModelDrummond("CO2")));
    }
    else if(params.activity_model == RTActivityModel::DebyeHuckel)
    {
        // Debye-Huckel full system
        aqueousphase.setActivityModel(chain(
                ActivityModelDebyeHuckelPHREEQC(),
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

    // Create a mineral phase
    MineralPhase mineralphase("Barite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(db, aqueousphase, mineralphase);

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

    double water_kg = 1.0;

    EquilibriumConditions conditions_ic(specs);
    conditions_ic.temperature(params.T, "celsius");
    conditions_ic.pressure(params.P, "bar");
    conditions_ic.pH(7.0); // H+ will be added as much as needed to achieve this pH

    // Define the initial condition (with formation water)
    ChemicalState state_ic(system);
    state_ic.setTemperature(params.T, "celsius");
    state_ic.setPressure(params.P, "bar");
    state_ic.setSpeciesMass("H2O(aq)", water_kg, "kg");
    state_ic.setSpeciesMass("SO4-2", 10 * water_kg, "ug");
    state_ic.setSpeciesMass("Ca+2", 995 * water_kg, "mg");
    state_ic.setSpeciesMass("Ba+2", 995 * water_kg, "mg");
    state_ic.setSpeciesMass("Sr+2", 105 * water_kg, "mg");
    state_ic.setSpeciesMass("Na+", 27250 * water_kg, "mg");
    state_ic.setSpeciesMass("K+", 1730 * water_kg, "mg");
    state_ic.setSpeciesMass("Mg+2", 110 * water_kg, "mg");
    state_ic.setSpeciesMass("Cl-", 45150 * water_kg, "mg");
    state_ic.setSpeciesMass("HCO3-", 1980 * water_kg, "mg");

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
    conditions_bc.pH(8.1); // H+ will be added as much as needed to achieve this pH

    // Define the initial condition (with seawater water)
    ChemicalState state_bc(system);
    state_bc.setTemperature(params.T, "celsius");
    state_bc.setPressure(params.P, "bar");
    state_bc.setSpeciesMass("H2O(aq)", water_kg, "kg");
    state_bc.setSpeciesMass("SO4-2", 2710 * water_kg, "mg");
    state_bc.setSpeciesMass("Ca+2", 411 * water_kg, "mg");
    state_bc.setSpeciesMass("Ba+2", 0.01 * water_kg, "mg");
    state_bc.setSpeciesMass("Sr+2", 8 * water_kg, "mg");
    state_bc.setSpeciesMass("Na+", 10760 * water_kg, "mg");
    state_bc.setSpeciesMass("K+", 399 * water_kg, "mg");
    state_bc.setSpeciesMass("Mg+2", 1290 * water_kg, "mg");
    state_bc.setSpeciesMass("Cl-", 19350 * water_kg, "mg");
    state_bc.setSpeciesMass("HCO3-", 142 * water_kg, "mg");

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
    output.add("speciesMolality(H+)");
    output.add("speciesMolality(Cl-)");
    output.add("speciesMolality(SO4-2)");
    output.add("speciesMolality(Ba+2)");
    output.add("speciesMolality(Ca+2)");
    output.add("speciesMolality(Sr+2)");
    output.add("speciesMolality(Na+)");
    output.add("speciesAmount(Barite)");
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
    output.add("elementmolality(Sr)");
    output.add("phaseAmount(Barite)");
    output.add("phaseMass(Barite)");
    output.add("phaseVolume(Barite)");
    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(REACTIVE_TRANSPORT_STEPS)

    // Reactive transport simulations in the cycle
    while (step <  params.nsteps)
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