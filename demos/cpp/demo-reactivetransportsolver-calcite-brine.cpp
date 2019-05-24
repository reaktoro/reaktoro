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
#include <sstream>      // for using stringstream
#include <iomanip>      // for setprecition
#include <tuple>

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
struct Params{

    // Discretisation params
    int ncells; // the number of cells in the spacial discretization
    int nsteps; // the number of steps in the reactive transport simulation
    double xl; // the x-coordinates of the left boundaries
    double xr; // the x-coordinates of the right boundaries
    double dx; // the space step (in units of m)
    double dt; // the time step (in units of s)

    // Physical params
    double D; // the diffusion coefficient (in units of m2/s)
    double v; // the Darcy velocity (in units of m/s)
    double T; // the temperature (in units of degC)
    double P; // the pressure (in units of bar)

    // Solver params
    bool is_smart_solver;
    double smart_reltol;
    double smart_abstol;

};

/// Forward declaration
auto mkdir(const std::string & folder) -> bool;
auto outputConsole(const Params & params) -> void;
auto makeResultsFolder(const Params & params) -> std::string;
auto runReactiveTransport(const Params & params) -> void;

int main()
{

    // Step 1: Initialise auxiliary time-related constants
    int second(1);
    int minute(60);
    int hour(60 * minute);
    int day(24 * hour);
    int year(365 * day);

    // Step 2: Define parameters for the reactive transport simulation
    Params params;

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 1.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    params.nsteps = 9600; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 2 * minute; // the time step (in units of s)

    // Define physical and chemical parameters
    params.D = 1.0e-9; // the diffusion coefficient (in units of m2/s)
    params.v = 1.0 / day; // the Darcy velocity (in units of m/s)
    params.T = 60.0;                     // the temperature (in units of degC)
    params.P = 100;                      // the pressure (in units of bar)

    // Define parameters of the equilirium solvers
    params.smart_reltol = 1e-1;
    params.smart_abstol = 1e-10;

    // Output
    outputConsole(params);

    // Execute reactive transport with different solvers
    params.is_smart_solver = 1; runReactiveTransport(params);
    params.is_smart_solver = 0; runReactiveTransport(params);

    return 0;
}
auto runReactiveTransport(const Params & params) -> void
{

    // Step **: Create the results folder
    auto folder = makeResultsFolder(params);

    // Step **: Define chemical equilibrium options
    EquilibriumOptions options;
    options.smart.reltol = params.smart_reltol;
    options.smart.abstol = params.smart_abstol;

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor;
    // Default chemical model (HKF extended Debye-Hückel model)
    editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO2(aq) CO3--");
    // Create aqueous phase with all possible elements
    // Set a chemical model of the phase with the Pitzer equation of state
    // With an exception for the CO2, for which Drummond model is set
    editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C")
            .setChemicalModelPitzerHMW()
            .setActivityModelDrummondCO2();
    editor.addMineralPhase("Quartz");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Dolomite");

    // Step **: Create the ChemicalSystem object using the configured editor
    ChemicalSystem system(editor);
    //if (params.is_smart_solver) std::cout << "system = \n" << system << std:: endl;

    // Step **: Define the initial condition (IC) of the reactive transport modeling problem
    EquilibriumProblem problem_ic(system);
    problem_ic.setTemperature(params.T, "celsius");
    problem_ic.setPressure(params.P, "bar");
    problem_ic.add("H2O",   1.0, "kg");
    problem_ic.add("NaCl",  0.7, "mol");
    problem_ic.add("CaCO3", 10,  "mol");
    problem_ic.add("SiO2",  10,  "mol");

    // Step **: Define the boundary condition (BC)  of the reactive transport modeling problem
    EquilibriumProblem problem_bc(system);
    problem_bc.setTemperature(params.T, "celsius");
    problem_bc.setPressure(params.P, "bar");
    problem_bc.add("H2O",   1.00, "kg");
    problem_bc.add("NaCl",  0.90, "mol");
    problem_bc.add("MgCl2", 0.05, "mol");
    problem_bc.add("CaCl2", 0.01, "mol");
    problem_bc.add("CO2",   0.75, "mol");

    // Step **: Calculate the equilibrium states for the IC and BC
    ChemicalState state_ic = equilibrate(problem_ic);
    ChemicalState state_bc = equilibrate(problem_bc);

    // Step **: Scale the volumes of the phases in the initial condition
    state_ic.scalePhaseVolume("Aqueous", 0.1, "m3");    //
    state_ic.scalePhaseVolume("Quartz", 0.882, "m3");
    state_ic.scalePhaseVolume("Calcite", 0.018, "m3");

    // Step **: Scale the boundary condition state
    state_bc.scaleVolume(1.0, "m3");

    // Step **: Create the mesh for the column
    Mesh mesh(params.ncells, params.xl, params.xr);

    // Step **: Create a chemical field object with every cell having state given by state_ic
    ChemicalField field(mesh.numCells(), state_ic);

    // Step **: Define the reactive transport modeling
    ReactiveTransportSolver rtsolver(system, params.is_smart_solver);
    rtsolver.setMesh(mesh);
    rtsolver.setVelocity(params.v);
    rtsolver.setDiffusionCoeff(params.D);
    rtsolver.setBoundaryState(state_bc);
    rtsolver.setTimeStep(params.dt);
    rtsolver.initialize();
    rtsolver.setEquilibriumOptions(options);

    // Step **: Define the quantities that should be output for every cell, every time step
    ChemicalOutput output(rtsolver.output());
    output.add("pH");
    output.add("speciesMolality(H+)");
    output.add("speciesMolality(Ca++)");
    output.add("speciesMolality(Mg++)");
    output.add("speciesMolality(HCO3-)");
    output.add("speciesMolality(CO2(aq))");
    output.add("phaseVolume(Calcite)");
    output.add("phaseVolume(Dolomite)");
    output.filename(folder + "/" + "test.txt");

    // Step **: Create two profilers on the step level and one general profiler
    Profiler rt_profiler(rtsolver.profile(Profiling::RT));
    Profiler eq_profiler(rtsolver.profile(Profiling::EQ));

    // Step **: Create a profiler on a cell level
    EquilibriumProfiler eq_cell_profiler(rtsolver.cellprofile(Profiling::EQ_CW));

    // Step **: Create status tracker
    SolverStatus tracker(rtsolver.trackStatus(folder, "status-tracker"));

    // Step **: Set initial time and counter of steps in time
    double t(0.0);
    int step(0);

    // Reactive transport simulations in the cycle
    while (step <= params.nsteps){

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step_tracked(field);

        // Increment time step and number of time steps
        t += params.dt;
        step += 1;
    }

    // Step **: Output the total time of the simulations
    rtsolver.outputProfiling();

    // Output profiling results
    rtsolver.outputProfiling(folder + "/profiling");
}

/// Make directory for Windows and Linux
auto mkdir(const std::string & folder) -> bool
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


/// Create results file with parameters of the test
auto makeResultsFolder(const Params & params) -> std::string
{
    struct stat status = {0};               // structure to get the file status

    std::ostringstream reltol_stream, abstol_stream, dt_stream;
    dt_stream << params.dt;
    reltol_stream << std::scientific << std::setprecision(1) << params.smart_reltol;
    abstol_stream << std::scientific << std::setprecision(1) <<  params.smart_abstol;

    std::string test_tag = "-dt-" + dt_stream.str() +
                           "-ncells-" + std::to_string(params.ncells) +
                           "-nsteps-" + std::to_string(params.nsteps) +
                           "-reltol-" + reltol_stream.str() +
                           "-abstol-" + abstol_stream.str() +
                           (params.is_smart_solver == true ? "-smart" : "-reference");      // name of the folder with results
    std::string folder = "../results" + test_tag;
    if (stat(folder.c_str(), &status) == -1) mkdir(folder.c_str());

    std::cout << "\nsolver                         : " << (params.is_smart_solver == true ? "smart" : "conventional") << std::endl;

    return folder;
}

auto outputConsole(const Params & params) -> void {

    // Log the parameters in the console
    std::cout << "dt      : " << params.dt << std::endl;
    std::cout << "ncells  : " << params.ncells << std::endl;
    std::cout << "nsteps  : " << params.nsteps << std::endl;
    std::cout << "D       : " << params.D << std::endl;
    std::cout << "v       : " << params.v << std::endl;
    std::cout << "CFD     : " << params.v * params.dt / params.dx << std::endl;
    std::cout << "T       : " << params.T << std::endl;
    std::cout << "P       : " << params.P << std::endl;
    std::cout << "abstol  : " << params.smart_abstol << std::endl;
    std::cout << "reltol  : " << params.smart_reltol << std::endl;

}

