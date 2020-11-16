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

// C++ includes
#include <sstream>      // for using stringstream
#include <iomanip>      // for setprecition

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

struct Params
{
    // Discretization params
    int ncells = 0; // the number of cells in the spacial discretization
    int nsteps = 0; // the number of steps in the reactive transport simulation
    double xl = 0; // the x-coordinates of the left boundaries
    double xr = 0; // the x-coordinates of the right boundaries
    double dx = 0; // the space step (in units of m)
    double dt = 0; // the time step (in units of s)

    // Physical params
    double D = 0; // the diffusion coefficient (in units of m2/s)
    double v = 0; // the Darcy velocity (in units of m/s)
    double T = 0; // the temperature (in units of degC)
    double P = 0; // the pressure (in units of bar)

    // Solver params
    bool use_smart_equilibrium_solver = false;
    double smart_equilibrium_reltol = 0;
    double amount_fraction_cutoff = 0;
    double mole_fraction_cutoff = 0;

    std::string activity_model = "";
};

struct Results
{
    /// Total CPU time (in s) required by smart equilibrium scheme
    double smart_total = 0.0;

    /// Total CPU time (in s) excluding the costs for the search of the closest reference states.
    double smart_total_ideal_search = 0.0;

    /// Total CPU time (in s) required by smart equilibrium scheme
    /// excluding the costs for the search and storage of the closest reference states.
    double smart_total_ideal_search_store = 0.0;

    /// Total CPU time (in s) required by conventional equilibrium scheme
    double conventional_total = 0.0;

    /// The total time taken to perform all time steps using conventional equilibrium algorithm
    double time_reactive_transport_conventional = 0.0;

    /// The total time taken to perform all time steps using smart equilibrium algorithm
    double time_reactive_transport_smart = 0.0;

    /// The accumulated timing information of all equilibrium calculations.
    EquilibriumTiming equilibrium_timing = {};

    /// The accumulated timing information of all smart equilibrium calculations.
    SmartEquilibriumTiming smart_equilibrium_timing = {};

    // Rate of the smart equilibrium estimation w.r.t to the total chemical equilibrium calculation
    double smart_equilibrium_acceptance_rate = 0.0;
};

/// Forward declaration
auto mkdir(const std::string& folder) -> bool;
auto outputConsole(const Params& params) -> void;
auto makeResultsFolder(const Params& params) -> std::string;
auto runReactiveTransport(const Params& params, Results& results) -> void;

int main()
{
    Time start = time();

    // Step 1: Initialise auxiliary time-related constants
    int minute = 60;
    int hour = 60 * minute;
    int day = 24 * hour;

    // Step 2: Define parameters for the reactive transport simulation
    Params params;

    // Define discretization parameters
    params.xl = 0.0; // the x-coordinates of the left boundaries
    params.xr = 100.0; // the x-coordinates of the right boundaries
    params.ncells = 100; // the number of cells in the spacial discretization
    params.nsteps = 5000; // the number of steps in the reactive transport simulation
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 0.1*day; // the time step (in units of s)



    // Define physical and chemical parameters
    params.D = 0.0;     // the diffusion coefficient (in units of m2/s)
    params.v = 3e-5; // the Darcy velocity (in units of m/s)
    params.T = 25.0;                     // the temperature (in units of degC)
    params.P = 1.01325;                      // the pressure (in units of bar)

    // Define parameters of the equilibrium solvers
    params.smart_equilibrium_reltol = 0.01;

    // Define the activity model for the aqueous species
    //params.activity_model = "dk-full";
    params.activity_model = "pitzer-full";
    //params.activity_model = "dk";

    // Define equilibrium solver cutoff tolerances
    params.amount_fraction_cutoff = 1e-14;
    params.mole_fraction_cutoff = 1e-14;

    // Output
    outputConsole(params);

    // Results
    Results results;

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
                                             - results.smart_equilibrium_timing.learning_storage;

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
    std::cout << "total time                          : " << elapsed(start) << std::endl;

    return 0;
}
auto runReactiveTransport(const Params& params, Results& results) -> void
{
    // Step **: Create the results folder
    auto folder = makeResultsFolder(params);

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
    if(params.activity_model == "hkf"){
        // HKF full system
        editor.addAqueousPhase(selected_species);
    }
    else if(params.activity_model == "pitzer"){
        // Pitzer selected species system
        editor.addAqueousPhase(selected_species)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "pitzer-full"){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    }
    else if(params.activity_model == "dk"){
        // Debye-Huckel selected species system
        editor.addAqueousPhase(selected_species)
                .setChemicalModelDebyeHuckel(dhModel);
    }
    else if(params.activity_model == "dk-full"){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel(dhModel);
    }
    editor.addMineralPhase("Siderite");
    editor.addMineralPhase("Pyrite");
    editor.addMineralPhase("Hematite");

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
    problem_ic.add("Siderite", 0.0, "mol");
    problem_ic.add("Pyrite", 0.0, "mol");
    problem_ic.add("Hematite", 0.5, "mol");
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
    problem_bc.add("Siderite", 0.0, "mol");
    problem_ic.add("Pyrite", 0.0, "mol");
    problem_ic.add("Hematite", 0.0, "mol");
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
    output.add("speciesAmount(Fe++)");
    output.add("speciesAmount(Siderite)");
    output.add("speciesAmount(Pyrite)");
    output.add("speciesAmount(Hematite)");
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
        // Print simulation progress
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }

    if(params.use_smart_equilibrium_solver)
        rtsolver.outputClusterInfo();

    if(params.use_smart_equilibrium_solver)
        results.time_reactive_transport_smart = toc(REACTIVE_TRANSPORT_STEPS);
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

/// Make directory for Windows and Linux
auto mkdir(const std::string& folder) -> bool
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
auto makeResultsFolder(const Params& params) -> std::string
{
    struct stat status = {0};               // structure to get the file status

    std::ostringstream reltol_stream, dt_stream;
    dt_stream << params.dt;
    reltol_stream << std::scientific << std::setprecision(1) << params.smart_equilibrium_reltol;

    std::string test_tag = "-dt-" + dt_stream.str() +
                           "-ncells-" + std::to_string(params.ncells) +
                           "-nsteps-" + std::to_string(params.nsteps) +
                           "-" + params.activity_model + "-reference";

    std::string smart_test_tag = "-dt-" + dt_stream.str() +
                                 "-ncells-" + std::to_string(params.ncells) +
                                 "-nsteps-" + std::to_string(params.nsteps) +
                                 "-reltol-" + reltol_stream.str() +
                                 "-" + params.activity_model +
                                 "-smart";

    std::string folder = "results-scavenging-with-hematite";
    folder = (params.use_smart_equilibrium_solver) ?
             folder + smart_test_tag :
             folder + test_tag;

    if (stat(folder.c_str(), &status) == -1) mkdir(folder);

    std::cout << "\nsolver                         : " << (params.use_smart_equilibrium_solver ? "smart" : "conventional") << std::endl;

    return folder;
}

auto outputConsole(const Params& params) -> void {

    // Log the parameters in the console
    std::cout << "dt      : " << params.dt << std::endl;
    std::cout << "ncells  : " << params.ncells << std::endl;
    std::cout << "nsteps  : " << params.nsteps << std::endl;
    std::cout << "D       : " << params.D << std::endl;
    std::cout << "v       : " << params.v << std::endl;
    std::cout << "CFD     : " << params.v * params.dt / params.dx << std::endl;
    std::cout << "T       : " << params.T << std::endl;
    std::cout << "P       : " << params.P << std::endl;
    std::cout << "eqreltol       : " << params.smart_equilibrium_reltol << std::endl;
    std::cout << "activity model : " << params.activity_model << std::endl;

}

