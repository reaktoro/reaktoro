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

// C++ includes
#include <algorithm>    // for using transform

#if defined _WIN32      // for creating a new folder
#include <windows.h>
#ifdef __MINGW32__
#include <sys/stat.h>
#endif
#else
#include <sys/stat.h>
#endif

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

using namespace Reaktoro;

/// Make directory for Windows and Linux
bool mkdir(const std::string & folder)
{
#if defined _WIN32
    // Replace slash by backslash
    std::transform(begin(folder), end(folder), begin(folder),
                   [](char ch) { return ch == '/' ? '\\' : ch; });
    return 0 != CreateDirectory(folder.c_str(), NULL);
#else
    // Create the directory with Read + Write + Execute rights for user, group, and others
    return ::mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif
}

int main()
{
    // Step 1: Initialise auxiliary time-related constants
    int second(1), minute(60);
    int hour(60 * minute), day(24 * hour), year(365 * day);

    // Step 2: Define parameters for the reactive transport simulation
    double xl(0.0), xr(100.0);          // the x-coordinates of the left and right boundaries
    int nsteps(100);                    // the number of steps in the reactive
    // transport simulation
    int ncells(100);                    // the number of cells in the spacial discretization
    double D(1.0e-9);                   // the diffusion coefficient (in units of m2/s)
    double v(1.0 / day);                // the fluid pore velocity (in units of m/s)
    double dt(0.5 * day);               // the time step (in units of s)
    double T(60.0);                     // the temperature (in units of degC)
    double P(100);                      // the pressure (in units of bar)
    double is_smrtsolv(true);          // the parameter that defines whether classic or smart EquilibriumSolver must be used

    // Step 3: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor;
    editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCl2 MgCl2 CO2");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Quartz");
    editor.addMineralPhase("Dolomite");

    // Step 4: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);

    // Step 5: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    problem_ic.setTemperature(T, "celsius");
    problem_ic.setPressure(P, "bar");
    problem_ic.add("H2O",   1.0, "kg");
    problem_ic.add("NaCl",  0.7, "mol");
    problem_ic.add("CaCO3", 10,  "mol");
    problem_ic.add("Si02",  10,  "mol");

    // Step 6: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumProblem problem_bc(system);
    problem_bc.setTemperature(T, "celsius");
    problem_bc.setPressure(P, "bar");
    problem_bc.add("H2O",   1.00, "kg");
    problem_bc.add("NaCl",  0.90, "mol");
    problem_bc.add("MgCl2", 0.05, "mol");
    problem_bc.add("CaCl2", 0.01, "mol");
    problem_bc.add("CO2",   0.75, "mol");

    // Step 7: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);


    // Step 8: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");
    state_ic.scalePhaseVolume("Quartz", 0.882, "m3");
    state_ic.scalePhaseVolume("Calcite", 0.018, "m3");

    // Step 9: Scale the boundary condition state
    state_bc.scaleVolume(1.0);

    // Step 10: Create the mesh for the column
    Mesh mesh(ncells, xl, xr);

    // Step 11: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step 12: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(system, is_smrtsolv);
    rtsolver.setMesh(mesh);
    rtsolver.setVelocity(v);
    rtsolver.setDiffusionCoeff(D);
    rtsolver.setBoundaryState(state_bc);
    rtsolver.setTimeStep(dt);
    rtsolver.initialize();

    // Step 13: Define the quantities that should be output for every cell, every time step
    ChemicalOutput output(rtsolver.output());
    output.add("pH");
    output.add("speciesMolality(H+)");
    output.add("speciesMolality(Ca++)");
    output.add("speciesMolality(Mg++)");
    output.add("speciesMolality(HCO3-)");
    output.add("speciesMolality(CO2(aq))");
    output.add("phaseVolume(Calcite)");
    output.add("phaseVolume(Dolomite)");

    // Step **: Create the results folder
    struct stat status = {0};               // structure to get the file status
    std::string test_tag = "-ncells-" + std::to_string(ncells) + "-nsteps-" + std::to_string(nsteps);      // name of the folder with results
    std::string folder = "../results" + test_tag;
    if (stat(folder.c_str(), &status) == -1) mkdir(folder.c_str()); // create the
    output.filename(folder + "/" + "test.txt");

    // Step **: Create two profilers
    Profiler rt_profiler(rtsolver.profile(Profiling::RT));
    Profiler eq_profiler(rtsolver.profile(Profiling::EQ));
    Profiler total_profiler(Profiling::Total);

    // Step 15: Set initial time and counter of steps in time
    double t(0.0);
    int step(0);

    // Step **: Start profiling total simulation
    total_profiler.startProfiling();

    // Reactive transport simulations in the cycle
    while (step <= nsteps){
        // Print the progress of the simulation
        std::cout << "Progress: " << step << " / " << nsteps << " "  << t/minute << " min" << std::endl;

        // Perform one reactive transport time step
        rtsolver.step(field, is_smrtsolv);

        // Increment time step and number of time steps
        t += dt;
        step += 1;
    }
    // Step **: End profiling total simulation
    total_profiler.endProfiling();

    // Step **: Output the total time of the simulations
    std::cout << "Elapsed time: ";
    total_profiler.consoleOutput();
    std::cout << std::endl;

    // Output profiling results
    rtsolver.outputProfiling(folder + "/profiling");
}
