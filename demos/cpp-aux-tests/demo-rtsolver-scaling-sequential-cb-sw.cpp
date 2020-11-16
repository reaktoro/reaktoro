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
    // Discretisation params
    int ncells = 0; // the number of cells in the spacial discretization
    double xl = 0; // the x-coordinates of the left boundaries
    double xr = 0; // the x-coordinates of the right boundaries
    double dx = 0; // the space step (in units of m)
    double dt = 0; // the time step (in units of s)

    int nsteps_cb = 0;  // the number of steps in the reactive transport simulation of the first injection phase
    int nsteps_sw = 0;  // the number of steps in the reactive transport simulation of the second injection phase
    int nsteps = 0;     // the total number of steps in the reactive transport simulation

    // Physical params
    double D = 0; // the diffusion coefficient (in units of m2/s)
    double v = 0; // the Darcy velocity (in units of m/s)
    double T = 0; // the temperature (in units of degC)
    double P = 0; // the pressure (in units of bar)
    double water_kg = 1.0;  // amount of water used in the experiment

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
    params.xr = 25.0; // the x-coordinates of the right boundaries
    params.ncells = 243; // the number of cells in the spacial discretization
    params.dx = (params.xr - params.xl) / params.ncells; // the time step (in units of s)
    params.dt = 1*hour; // the time step (in units of s)

    params.nsteps_cb = 45;  // the number of steps in the reactive transport simulation of the first injection phase
    params.nsteps_sw = 855;  // the number of steps in the reactive transport simulation of the second injection phase
    params.nsteps = params.nsteps_cb + params.nsteps_sw;     // the total number of steps in the reactive transport simulation

    // Define physical and chemical parameters
    params.D = 0.0;             // the diffusion coefficient (in units of m2/s)
    params.v = 0.8e-5;          // the Darcy velocity (in units of m/s)
    params.T = 60.0;            // the temperature (in units of degC)
    params.P = 200 * 1.01325;   // the pressure (in units of bar)

    // Define parameters of the equilibrium solvers
    params.smart_equilibrium_reltol = 0.003;

    // Define the activity model for the aqueous species
    params.activity_model = "dk-full";
    //params.activity_model = "pitzer-full";
    //params.activity_model = "hkf-full";

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

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    // Step **: Construct the chemical system with its phases and species (using ChemicalEditor)
    ChemicalEditor editor(database);

    // Define the list of selected elements
    StringList selected_elements = "H Cl S O Ba Ca Sr Na K Mg C Si";

    // Depending on the activity model, define it using ChemicalEditor
    if(params.activity_model == "dk-full"){
        // Debye-Huckel full system
        editor.addAqueousPhaseWithElements(selected_elements)
                .setChemicalModelDebyeHuckel(dhModel);
    }
    else if(params.activity_model == "pitzer-full"){
        // Debye-Huckel full system
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
    problem_ic.add("H2O", params.water_kg, "kg");
    problem_ic.add("SO4", 10 * params.water_kg, "ug");
    problem_ic.add("Ca", 995 * params.water_kg, "mg");
    problem_ic.add("Ba", 995 * params.water_kg, "mg");
    problem_ic.add("Sr", 105 * params.water_kg, "mg");
    problem_ic.add("Na", 27250 * params.water_kg, "mg");
    problem_ic.add("K", 1730 * params.water_kg, "mg");
    problem_ic.add("Mg", 110 * params.water_kg, "mg");
    problem_ic.add("Cl", 45150 * params.water_kg, "mg");
    problem_ic.add("HCO3", 1980 * params.water_kg, "mg");
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
    problem_bc_cb.add("H2O", params.water_kg, "kg");
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
    problem_bc_sw.add("H2O", params.water_kg, "kg");
    problem_bc_sw.add("SO4--", 2710 * params.water_kg, "mg");
    problem_bc_sw.add("Ca++", 411 * params.water_kg, "mg");
    problem_bc_sw.add("Ba++", 0.01 * params.water_kg, "mg");
    problem_bc_sw.add("Sr++", 8 * params.water_kg, "mg");
    problem_bc_sw.add("Na+", 10760 * params.water_kg, "mg");
    problem_bc_sw.add("K+", 399 * params.water_kg, "mg");
    problem_bc_sw.add("Mg++", 1290 * params.water_kg, "mg");
    problem_bc_sw.add("Cl-", 19350 * params.water_kg, "mg");
    problem_bc_sw.add("HCO3-", 142 * params.water_kg, "mg");
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
    rtsolver.setBoundaryState(state_bc_cb);
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
    output.filename(folder + "/" + "test.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler;

    // Step **: Set initial time and counter of steps in time
    double t = 0.0;
    int step = 0;

    tic(REACTIVE_TRANSPORT_STEPS);

    // Reactive transport simulations in the cycle
    while (step <  params.nsteps_cb)
    {
        // // Define activity model depending on the parameter
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver.step(field);

        // Update the profiler after every call to step method
        profiler.update(rtsolver.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }

    if(params.use_smart_equilibrium_solver){
        std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
        std::cout << " CLUSTERS: " << std::endl;
        std::cout << "--------------------------------------------------------------------------------------------------" << std::endl;
        rtsolver.outputClusterInfo();
    }

    // Step **: Define new reactive transport instance for injecting seawater
    ReactiveTransportSolver rtsolver_sw(system);
    rtsolver_sw.setOptions(reactive_transport_options);
    rtsolver_sw.setMesh(mesh);
    rtsolver_sw.setVelocity(params.v);
    rtsolver_sw.setDiffusionCoeff(params.D);
    rtsolver_sw.setBoundaryState(state_bc_sw);
    rtsolver_sw.setTimeStep(params.dt);
    rtsolver_sw.initialize();

    // Step **: Define the quantities that should be output for every cell, every time step
    ChemicalOutput output_sw(rtsolver_sw.output());
    output_sw.add("pH");
    output_sw.add("speciesMolality(H+)");
    output_sw.add("speciesMolality(Cl-)");
    output_sw.add("speciesMolality(SO4--)");
    output_sw.add("speciesMolality(Ba++)");
    output_sw.add("speciesMolality(Ca++)");
    output_sw.add("speciesMolality(Sr++)");
    output_sw.add("speciesMolality(Na+)");
    output_sw.add("speciesMolality(Barite)");
    output_sw.add("elementmolality(Ba)");
    output_sw.add("elementmolality(C)");
    output_sw.add("elementmolality(Ca)");
    output_sw.add("elementmolality(Cl)");
    output_sw.add("elementmolality(H)");
    output_sw.add("elementmolality(K)");
    output_sw.add("elementmolality(Mg)");
    output_sw.add("elementmolality(Na)");
    output_sw.add("elementmolality(O)");
    output_sw.add("elementmolality(S)");
    output_sw.add("elementmolality(Si)");
    output_sw.add("elementmolality(Sr)");
    output_sw.add("elementmolality(Z)");
    output.add("phaseAmount(Barite)");
    output.add("phaseMass(Barite)");
    output.add("phaseVolume(Barite)");

    output_sw.filename(folder + "/" + "test-sw.txt");

    // Step **: Create RTProfiler to track the timing and results of reactive transport
    ReactiveTransportProfiler profiler_sw;

    t = 0.0;
    step = 0;

    while (step < params.nsteps_sw)
    {
        // // Define activity model depending on the parameter
        std::cout << "Step " << step << " of " << params.nsteps << std::endl;

        // Perform one reactive transport time step (with profiling of some parts of the transport simulations)
        rtsolver_sw.step(field);

        // Update the profiler after every call to step method
        profiler_sw.update(rtsolver_sw.result());

        // Increment time step and number of time steps
        t += params.dt;

        step += 1;
    }

    if(params.use_smart_equilibrium_solver){
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " CLUSTERS: " << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        rtsolver_sw.outputClusterInfo();
    }

    if(params.use_smart_equilibrium_solver) results.time_reactive_transport_smart = toc(REACTIVE_TRANSPORT_STEPS);
    else results.time_reactive_transport_conventional = toc(REACTIVE_TRANSPORT_STEPS);

    // Step **: Collect the analytics related to reactive transport performance
    auto analysis = profiler.analysis();
    auto rt_results = profiler.results();
    // Step **: Collect the analytics related to reactive transport performance
    auto analysis_sw = profiler_sw.analysis();
    auto rt_results_sw = profiler_sw.results();

    // Step **: Generate json output file with collected profiling data
    if(params.use_smart_equilibrium_solver)  JsonOutput(folder + "/" + "analysis-smart.json") << analysis;
    else    JsonOutput(folder + "/" + "analysis-conventional.json") << analysis;

    // Step **: Generate json output file with collected profiling data
    if(params.use_smart_equilibrium_solver)  JsonOutput(folder + "/" + "analysis-smart-sw.json") << analysis_sw;
    else    JsonOutput(folder + "/" + "analysis-conventional-sw.json") << analysis_sw;

    // Step **: Save equilibrium timing to compare the speedup of smart equilibrium solver versus conventional one
    if(params.use_smart_equilibrium_solver) {
        results.smart_equilibrium_timing = analysis.smart_equilibrium.timing;
        results.smart_equilibrium_acceptance_rate = analysis.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;

        std::cout << "smart equilibrium acceptance rate (cb) : " << results.smart_equilibrium_acceptance_rate << " / "
                  << (1 - results.smart_equilibrium_acceptance_rate) * params.ncells * params.nsteps_cb
                  << " fully evaluated GEMS out of " << params.ncells * params.nsteps_cb  << std::endl;

    }
    else results.equilibrium_timing = analysis.equilibrium.timing;

    // Step **: Save equilibrium timing to compare the speedup of smart equilibrium solver versus conventional one
    if(params.use_smart_equilibrium_solver) {
        results.smart_equilibrium_timing += analysis_sw.smart_equilibrium.timing;
        results.smart_equilibrium_acceptance_rate = analysis_sw.smart_equilibrium.smart_equilibrium_estimate_acceptance_rate;

        std::cout << "smart equilibrium acceptance rate (sw) : " << results.smart_equilibrium_acceptance_rate << " / "
                  << (1 - results.smart_equilibrium_acceptance_rate) * params.ncells * params.nsteps_sw
                  << " fully evaluated GEMS out of " << params.ncells * params.nsteps_sw  << std::endl;

    }
    else results.equilibrium_timing = analysis_sw.equilibrium.timing;
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

    std::string folder = "results-scaling-sequential-cb-sw";
    folder = (params.use_smart_equilibrium_solver) ?
             folder + smart_test_tag :
             folder + test_tag;

    if (stat(folder.c_str(), &status) == -1) mkdir(folder.c_str());

    std::cout << "\nsolver                         : " << (params.use_smart_equilibrium_solver == true ? "smart" : "conventional") << std::endl;

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

