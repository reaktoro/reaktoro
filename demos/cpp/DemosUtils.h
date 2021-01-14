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

#pragma once

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

/// The options for the activity model of aqueous species
enum ActivityModel
{
    HKF,

    DebyeHuckel,

    Pitzer,

    HKFSelectedSpecies,

    DebyeHuckelSelectedSpecies,

    PitzerSelectedSpecies,

};

// Forward declaration
auto getSmartMethodTag(enum SmartEquilibriumStrategy method) -> std::string;
auto getSmartMethodTag(enum SmartKineticStrategy method) -> std::string;
auto getActivityModelTag(enum ActivityModel activity_model) -> std::string;
auto getGibbsHessianTag(enum GibbsHessian hessian) -> std::string;
auto mkdir(std::string& folder) -> bool;

struct ReactiveTransportParams
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

    // ODML params
    bool use_smart_equilibrium_solver = false;
    bool use_smart_kinetics_solver = false;

    double smart_equilibrium_reltol = 0;
    double amount_fraction_cutoff = 0;
    double mole_fraction_cutoff = 0;

    double smart_kinetic_tol = 0.0;
    double smart_kinetic_reltol = 0.0;
    double smart_kinetic_abstol = 0.0;

    ActivityModel activity_model = ActivityModel::HKF;

    GibbsHessian hessian = GibbsHessian::Exact;

    SmartEquilibriumStrategy method = SmartEquilibriumStrategy::Clustering;
    SmartKineticStrategy smart_kin_method = SmartKineticStrategy::Clustering;

    // Outputting params
    double output_results = true;

    /// Output test parameters to the console for the reactive transport with equilibrium only
    auto outputConsole() -> void
    {
        // Log the parameters in the console
        std::cout << "dt      : " << dt << std::endl;
        std::cout << "ncells  : " << ncells << std::endl;
        std::cout << "nsteps  : " << nsteps << std::endl;
        std::cout << "D       : " << D << std::endl;
        std::cout << "v       : " << v << std::endl;
        std::cout << "CFD     : " << v * dt / dx << std::endl;
        std::cout << "T       : " << T << std::endl;
        std::cout << "P       : " << P << std::endl;
        std::cout << "activity model  : " << getActivityModelTag(activity_model) << std::endl;
        std::cout << "hessian         : " << getGibbsHessianTag(hessian) << std::endl;

    }

    /// Output test parameters to the console for the reactive transport with equilibrium only
    auto outputEquilibriumMethod() -> void
    {
        std::cout << "smart eq used   : " << use_smart_equilibrium_solver << std::endl;
        std::cout << "smart eq method : " << getSmartMethodTag(method) << std::endl;
        std::cout << "eqreltol        : " << smart_equilibrium_reltol << std::endl;
    }

    /// Output test parameters to the console for the reactive transport with kinetics only
    auto outputConsoleKineticMethod() -> void
    {
        // Log the parameters in the console
        outputEquilibriumMethod();
        std::cout << "smart kin used   : " << use_smart_kinetics_solver << std::endl;
        std::cout << "smart kin method : " << getSmartMethodTag(smart_kin_method) << std::endl;
        std::cout << "kinetics tol     : " << smart_kinetic_tol << std::endl;
        std::cout << "kinetics reltol  : " << smart_kinetic_reltol << std::endl;
        std::cout << "kinetics abstol  : " << smart_kinetic_abstol << std::endl;
    }

    /// Create results file with parameters of the test for the reactive transport with equilibrium only
    auto makeResultsFolder(const std::string& demo_tag) -> std::string
    {
        struct stat status = {0};               // structure to get the file status

        std::ostringstream reltol_stream, dt_stream;
        dt_stream << dt;
        reltol_stream << std::scientific << std::setprecision(1) << smart_equilibrium_reltol;

        std::string test_tag = "-dt-" + dt_stream.str() +
                               "-ncells-" + std::to_string(ncells) +
                               "-nsteps-" + std::to_string(nsteps) +
                               "-" + getActivityModelTag(activity_model) + "-conventional";

        std::string smart_test_tag = "-" + getSmartMethodTag(method) +
                                     "-dt-" + dt_stream.str() +
                                     "-ncells-" + std::to_string(ncells) +
                                     "-nsteps-" + std::to_string(nsteps) +
                                     "-reltol-" + reltol_stream.str() +
                                     "-" + getActivityModelTag(activity_model) + "-smart";

        std::string folder = "rt-" + demo_tag;
        folder = (use_smart_equilibrium_solver) ?
                 folder + smart_test_tag :
                 folder + test_tag;

        if (stat(folder.c_str(), &status) == -1) mkdir(folder);

        std::cout << "\nsolver                         : " << (use_smart_equilibrium_solver ? "smart" : "conventional") << std::endl;

        return folder;
    }

    /// Create results file with parameters of the test for the reactive transport with kinetics
    auto makeResultsFolderKinetics(const std::string& demo_tag) -> std::string
    {
        struct stat status = {0};               // structure to get the file status

        std::ostringstream eqreltol_stream,
                eqabstol_stream,
                eqcutoff_stream,
                eqtol_stream,
                dt_stream,
                kinreltol_stream,
                kinabstol_stream,
                kintol_stream;
        dt_stream << dt;
        eqreltol_stream << std::scientific << std::setprecision(1) << smart_equilibrium_reltol;
        kinreltol_stream << std::scientific << std::setprecision(1) << smart_kinetic_reltol;
        kinabstol_stream << std::scientific << std::setprecision(1) << smart_kinetic_abstol;
        kintol_stream << std::scientific << std::setprecision(1) << smart_kinetic_tol;
        std::string test_tag = "-dt-" + dt_stream.str() +
                               "-ncells-" + std::to_string(ncells) +
                               "-nsteps-" + std::to_string(nsteps) +
                               "-" + getActivityModelTag(activity_model) +
                               (use_smart_kinetics_solver ? "-smart-kin" : "-conv-kin") +
                               (use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq");      // name of the folder with results
        std::string smart_test_tag = "-" + getSmartMethodTag(method) +
                                     "-dt-" + dt_stream.str() +
                                     "-ncells-" + std::to_string(ncells) +
                                     "-nsteps-" + std::to_string(nsteps) +
                                     "-eqrel-" + eqreltol_stream.str() +
                                     "-kinrel-" + kinreltol_stream.str() +
                                     "-kinabs-" + kinabstol_stream.str() +
                                     "-kintol-" + kintol_stream.str() +
                                     "-" + getActivityModelTag(activity_model) +
                                     (use_smart_kinetics_solver ? "-smart-kin" : "-conv-kin") +
                                     (use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq");      // name of the folder with results

        std::string folder = "rt-kinetics-" + demo_tag;
        folder = (use_smart_kinetics_solver || use_smart_equilibrium_solver) ?
                 demo_tag + smart_test_tag :
                 demo_tag + test_tag;
        if (stat(folder.c_str(), &status) == -1) mkdir(folder);

        std::cout << "*********************************************************************" << std::endl;
        std::cout << "*********************************************************************" << std::endl;
        std::cout << "\nsolver                         : "
                  << (use_smart_kinetics_solver ? "smart_kin & " : "conv_kin & ")
                  << (use_smart_equilibrium_solver ? "smart_eq" : "conv_eq") << std::endl;

        return folder;
    }

};

struct ReactiveTransportResults
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

struct ReactiveTransportKineticsResults
{
    // Conventional kinetic and conventional equilibrium schemes' times
    // *********************************************************************************//

    /// Total CPU time (in s) required by conventional kinetic and equilibrium schemes.
    double conv_kin_conv_eq_total = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and equilibrium schemes
    /// excluding the costs for the chemical properties evaluation.
    double conv_kin_conv_eq_total_ideal_properties = 0.0;

    /// Total CPU time (in s) required for equilibrium in the conventional kinetic using equilibrium schemes
    double conv_kin_conv_eq_total_equilibration = 0.0;

    // Smart kinetic and conventional equilibrium schemes' times
    // *********************************************************************************//

    /// Total CPU time (in s) required by smart kinetic and conventional equilibrium schemes.
    double smart_kin_conv_eq_total = 0.0;

    /// Total CPU time (in s) required by smart kinetic and conventional equilibrium schemes
    /// excluding the costs for the search of the closest reference states.
    double smart_kin_conv_eq_total_ideal_search = 0.0;

    /// Total CPU time (in s) required by smart kinetic and conventional equilibrium schemes
    /// excluding the costs for the search and storage of the closest reference states.
    double smart_kin_conv_eq_total_ideal_search_store = 0.0;

    /// Total CPU time (in s) required by smart kinetic and conventional equilibrium schemes
    /// excluding the costs for the search and storage of the closest reference states
    /// and the chemical properties evaluation.
    double smart_kin_conv_eq_total_ideal_search_store_properties = 0.0;

    // Rate of the smart kinetic estimation w.r.t to the total chemical kinetics calculation
    double smart_kin_conv_eq_acceptance_rate = 0.0;

    // Conventional kinetic and smart equilibrium schemes' times
    // *********************************************************************************//

    double conv_kin_smart_eq_total = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search of the closest reference states.
    double conv_kin_smart_eq_total_ideal_search = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search and store of the closest reference states.
    double conv_kin_smart_eq_total_ideal_search_store = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search and store of the closest reference states
    /// and the chemical properties evaluation.
    double conv_kin_smart_eq_total_ideal_search_store_properties = 0.0;


    /// Total CPU time (in s) required for smart equilibrium in the conventional kinetic using smart equilibrium schemes
    double conv_kin_smart_eq_total_smart_equilibration = 0.0;

    // Rate of the smart equilibrium estimation w.r.t to the total chemical kinetics calculation
    double conv_kin_smart_eq_equilibrium_acceptance_rate = 0.0;

    // Smart kinetic and smart equilibrium schemes' times
    // *********************************************************************************//

    double smart_kin_smart_eq_total = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search of the closest reference states.
    double smart_kin_smart_eq_total_ideal_search = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search and store of the closest reference states.
    double smart_kin_smart_eq_total_ideal_search_store = 0.0;

    /// Total CPU time (in s) required by conventional kinetic and smart equilibrium schemes
    /// excluding the costs for the search and store of the closest reference states
    /// and the chemical properties evaluation.
    double smart_kin_smart_eq_total_ideal_search_store_properties = 0.0;

    // Rate of the smart kinetics estimation w.r.t to the total chemical kinetics calculation
    double smart_kin_smart_eq_acceptance_rate = 0.0;

    // Rate of the smart equilibrium estimation w.r.t to the total chemical kinetics calculation
    double  smart_kin_smart_eq_equilibrium_acceptance_rate = 0.0;

    // Total times of reactive transport steps
    // *********************************************************************************//

    /// Total time taken to perform all time steps using conventional kinetics and conventional equilibrium algorithm
    double time_reactive_transport_conv_kin_conv_eq = 0.0;

    /// Total time taken to perform all time steps using smart kinetics and conventional equilibrium algorithm
    double time_reactive_transport_smart_kin_conv_eq = 0.0;

    /// Total time taken to perform all time steps using conventional kinetics and smart equilibrium algorithm
    double time_reactive_transport_conv_kin_smart_eq = 0.0;

    /// Total time taken to perform all time steps using smart kinetics and smart equilibrium algorithm
    double time_reactive_transport_smart_kin_smart_eq = 0.0;

    // Total times of reactive transport steps
    // *********************************************************************************//

    /// Accumulated timing information of all equilibrium calculations.
    EquilibriumTiming equilibrium_timing = {};

    /// The accumulated timing information of all smart equilibrium calculations.
    SmartEquilibriumTiming smart_equilibrium_timing = {};

    /// Accumulated timing information of all kinetic calculations.
    KineticTiming kinetic_timing = {};

    /// The accumulated timing information of all smart kinetic calculations.
    SmartKineticTiming smart_kinetic_timing = {};

    auto outputStatisticsSmartKineticsConventionalEquilibrium(const ReactiveTransportParams& params) -> void
    {
        smart_kin_conv_eq_total = smart_kinetic_timing.solve;
        smart_kin_conv_eq_total_ideal_search = smart_kinetic_timing.solve
                                                       - smart_kinetic_timing.estimate_search;
        smart_kin_conv_eq_total_ideal_search_store = smart_kinetic_timing.solve
                                                             - smart_kinetic_timing.estimate_search
                                                             - smart_kinetic_timing.learn_storage;
        smart_kin_conv_eq_total_ideal_search_store_properties = smart_kinetic_timing.solve
                                                                        - smart_kinetic_timing.learn_chemical_properties
                                                                        - smart_kinetic_timing.estimate_search
                                                                        - smart_kinetic_timing.learn_storage;

        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " - solve                : " << smart_kinetic_timing.solve << std::endl;
        std::cout << "   - learn              : " << smart_kinetic_timing.learn << " (" << smart_kinetic_timing.learn / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - integrate             : " << smart_kinetic_timing.learn_integration << " (" << smart_kinetic_timing.learn_integration / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - chemical properties : " << smart_kinetic_timing.learn_chemical_properties << " (" << smart_kinetic_timing.learn_chemical_properties / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - integrate_reaction_rates      : " << smart_kinetic_timing.learn_reaction_rates << " (" << smart_kinetic_timing.learn_reaction_rates / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - equilibration       : " << smart_kinetic_timing.learn_equilibration << " (" << smart_kinetic_timing.learn_equilibration / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - store                 : " << smart_kinetic_timing.learn_storage << " (" << smart_kinetic_timing.learn_storage / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "   - estimate           : " << smart_kinetic_timing.estimate << " (" << smart_kinetic_timing.estimate / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - search                : " << smart_kinetic_timing.estimate_search << " (" << smart_kinetic_timing.estimate_search / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - acceptance            : " << smart_kinetic_timing.estimate_error_control << " (" << smart_kinetic_timing.estimate_error_control / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "   - equilibrate           : " << smart_kinetic_timing.equilibrate << " (" << smart_kinetic_timing.equilibrate / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " acceptance rate      : " << smart_kin_conv_eq_acceptance_rate << " / " << (1 - smart_kin_conv_eq_acceptance_rate) * params.ncells * params.nsteps << " learnings out of " << params.ncells * params.nsteps  << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " - solve - search                 : " << smart_kin_conv_eq_total_ideal_search << std::endl;
        std::cout << " - solve - search - store         : " << smart_kin_conv_eq_total_ideal_search_store << std::endl;
        std::cout << " - solve - search - store - prop. : " << smart_kin_conv_eq_total_ideal_search_store_properties << std::endl;

    }

    auto outputStatisticsSmartKineticsSmartEquilibrium(const ReactiveTransportParams& params) -> void
    {
        smart_kin_smart_eq_total = smart_kinetic_timing.solve;
        smart_kin_smart_eq_total_ideal_search = smart_kinetic_timing.solve
                                                        - smart_kinetic_timing.estimate_search
                                                        - smart_equilibrium_timing.estimate_search;
        smart_kin_smart_eq_total_ideal_search_store = smart_kinetic_timing.solve
                                                              - smart_kinetic_timing.estimate_search
                                                              - smart_kinetic_timing.learn_storage
                                                              - smart_equilibrium_timing.estimate_search
                                                              - smart_equilibrium_timing.learn_storage;
        smart_kin_smart_eq_total_ideal_search_store_properties = smart_kinetic_timing.solve
                                                                         - smart_kinetic_timing.learn_chemical_properties
                                                                         - smart_kinetic_timing.estimate_search
                                                                         - smart_kinetic_timing.learn_storage
                                                                         - smart_equilibrium_timing.learn_chemical_properties
                                                                         - smart_equilibrium_timing.estimate_search
                                                                         - smart_equilibrium_timing.learn_storage;

        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " - solve                : " << smart_kinetic_timing.solve << std::endl;
        std::cout << "   - learn              : " << smart_kinetic_timing.learn << " (" << smart_kinetic_timing.learn / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - integrate             : " << smart_kinetic_timing.learn_integration << " (" << smart_kinetic_timing.learn_integration / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - chemical properties      : " << smart_kinetic_timing.learn_chemical_properties << " (" << smart_kinetic_timing.learn_chemical_properties / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - integrate_reaction_rates : " << smart_kinetic_timing.learn_reaction_rates << " (" << smart_kinetic_timing.learn_reaction_rates / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - equilibration            : " << smart_kinetic_timing.learn_equilibration << " (" << smart_kinetic_timing.learn_equilibration / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - store                 : " << smart_kinetic_timing.learn_storage << " (" << smart_kinetic_timing.learn_storage / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "   - estimate           : " << smart_kinetic_timing.estimate << " (" << smart_kinetic_timing.estimate / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - search                : " << smart_kinetic_timing.estimate_search << " (" << smart_kinetic_timing.estimate_search / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - acceptance            : " << smart_kinetic_timing.estimate_error_control << " (" << smart_kinetic_timing.estimate_error_control / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "   - equilibrate        : " << smart_kinetic_timing.equilibrate << " (" << smart_kinetic_timing.equilibrate / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - learning              : " << smart_equilibrium_timing.learn << " (" << smart_equilibrium_timing.learn / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - store                    : " << smart_equilibrium_timing.learn_storage << " (" << smart_equilibrium_timing.learn_storage / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - estimation            : " << smart_equilibrium_timing.estimate << " (" << smart_equilibrium_timing.estimate / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - search                   : " << smart_equilibrium_timing.estimate_search << " (" << smart_equilibrium_timing.estimate_search / smart_kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " smart kinetics acceptance rate      : " << smart_kin_smart_eq_acceptance_rate << " ( "
                  << smart_kin_smart_eq_acceptance_rate * 100 << " % ) / "
                  << (1 - smart_kin_smart_eq_acceptance_rate) * params.ncells * params.nsteps
                  << " learnings out of " << params.ncells * params.nsteps
                  << " ( " << (1 - smart_kin_smart_eq_acceptance_rate) * 100 << "% )" << std::endl;
        std::cout << " smart equilibrium acceptance rate   : " << smart_kin_smart_eq_equilibrium_acceptance_rate << " ( "
                  << smart_kin_smart_eq_equilibrium_acceptance_rate * 100 << " % ) / "
                  << (1 - smart_kin_smart_eq_equilibrium_acceptance_rate) * params.ncells * params.nsteps
                  << " learnings states out of " << params.ncells * params.nsteps
                  << " ( " << (1 - smart_kin_smart_eq_equilibrium_acceptance_rate) * 100 << "% )" << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " - solve - search                 : " << smart_kin_smart_eq_total_ideal_search << std::endl;
        std::cout << " - solve - search - store         : " << smart_kin_smart_eq_total_ideal_search_store << std::endl;
        std::cout << " - solve - search - store - prop. : " << smart_kin_smart_eq_total_ideal_search_store_properties << std::endl;
    }

    auto outputStatisticsConventionalKineticsConventionalEquilibrium(const ReactiveTransportParams& params) -> void
    {
        conv_kin_conv_eq_total = kinetic_timing.solve;
        conv_kin_conv_eq_total_ideal_properties = kinetic_timing.solve - kinetic_timing.integrate_chemical_properties;
        conv_kin_conv_eq_total_equilibration = kinetic_timing.integrate_equilibration;

        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " - solve                   : " << kinetic_timing.solve << std::endl;
        std::cout << "   - initialize            : " << kinetic_timing.initialize << " (" << kinetic_timing.initialize / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "   - integrate             : " << kinetic_timing.integrate << " (" << kinetic_timing.integrate / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - chemical properties      : " << kinetic_timing.integrate_chemical_properties << " (" << kinetic_timing.integrate_chemical_properties / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - integrate_reaction_rates : " << kinetic_timing.integrate_reaction_rates << " (" << kinetic_timing.integrate_reaction_rates / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - equilibration            : " << kinetic_timing.integrate_equilibration << " (" << kinetic_timing.integrate_equilibration / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "   - equilibrate           : " << kinetic_timing.equilibrate << " (" << kinetic_timing.equilibrate / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
    }

    auto outputStatisticsConventionalKineticsSmartEquilibrium(const ReactiveTransportParams& params) -> void
    {
        conv_kin_smart_eq_total = kinetic_timing.solve;
        conv_kin_smart_eq_total_ideal_search = kinetic_timing.solve - smart_equilibrium_timing.estimate_search;
        conv_kin_smart_eq_total_ideal_search_store = kinetic_timing.solve
                                                             - smart_equilibrium_timing.estimate_search
                                                             - smart_equilibrium_timing.learn_storage;
        conv_kin_smart_eq_total_ideal_search_store_properties = kinetic_timing.solve
                                                                        - smart_equilibrium_timing.estimate_search
                                                                        - smart_equilibrium_timing.learn_storage
                                                                        - kinetic_timing.integrate_chemical_properties;
        conv_kin_smart_eq_total_smart_equilibration = kinetic_timing.equilibrate;

        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " - solve                   : " << kinetic_timing.solve << std::endl;
        std::cout << "   - initialize            : " << kinetic_timing.initialize << " (" << kinetic_timing.initialize / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "   - integrate             : " << kinetic_timing.integrate << " (" << kinetic_timing.integrate / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - chemical properties : " << kinetic_timing.integrate_chemical_properties << " (" << kinetic_timing.integrate_chemical_properties / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - integrate_reaction_rates      : " << kinetic_timing.integrate_reaction_rates << " (" << kinetic_timing.integrate_reaction_rates / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - equilibration       : " << kinetic_timing.equilibrate << " (" << kinetic_timing.equilibrate / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "   - equilibrate             : " << kinetic_timing.equilibrate << " (" << kinetic_timing.equilibrate / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - learning          : " << smart_equilibrium_timing.learn << " (" << smart_equilibrium_timing.learn / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - store           : " << smart_equilibrium_timing.learn_storage << " (" << smart_equilibrium_timing.learn_storage / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "     - estimation        : " << smart_equilibrium_timing.estimate << " (" << smart_equilibrium_timing.estimate / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "       - search          : " << smart_equilibrium_timing.estimate_search << " (" << smart_equilibrium_timing.estimate_search / kinetic_timing.solve * 100 << " %)" << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " smart equilibrium acceptance rate   : " << conv_kin_smart_eq_equilibrium_acceptance_rate << " ( "
                  << conv_kin_smart_eq_equilibrium_acceptance_rate * 100 << " % ) / "
                  << (1 - conv_kin_smart_eq_equilibrium_acceptance_rate) * params.ncells * params.nsteps
                  << " learnings states out of " << params.ncells * params.nsteps
                  << " ( " << (1 - conv_kin_smart_eq_equilibrium_acceptance_rate) * 100 << "% )" << std::endl;
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << " - solve - search              : " << kinetic_timing.solve
                                                            - smart_equilibrium_timing.estimate_search << std::endl;
        std::cout << " - solve - search - store      : " << kinetic_timing.solve
                                                            - smart_equilibrium_timing.estimate_search
                                                            - smart_equilibrium_timing.learn_storage << std::endl;
        std::cout << " - solve - search - store - properties : " << kinetic_timing.solve
                                                                    - kinetic_timing.integrate_chemical_properties
                                                                    - smart_equilibrium_timing.learn_storage
                                                                    - smart_equilibrium_timing.estimate_search << std::endl;
    }


};

struct KineticPathParams
{
    // Discretization params
    double t0 = 0;      // starting time of simulations
    double tfinal = 0;  // final time of simulation

    double dt = 0.0;    // discretization step in time
    int n = 0;          // number of steps

    // Thermodynamic parameters
    double T = 0; // the temperature (in units of degC)
    double P = 0; // the pressure (in units of bar)

    ActivityModel activity_model = ActivityModel::HKF;

    /// -------------------------------------------------------------------------- ///
    /// Smart equilibrium calculations
    /// -------------------------------------------------------------------------- ///

    // ODML parameters for smart equilibrium calculations
    bool use_smart_equilibrium_solver = false;

    // Smart equilibrium method with clustering
    SmartEquilibriumStrategy smart_eq_method = SmartEquilibriumStrategy::Clustering;
    double smart_equilibrium_reltol = 1e-3;

    /// -------------------------------------------------------------------------- ///
    /// Smart kinetics calculations
    /// -------------------------------------------------------------------------- ///

    // ODML parameters for smart kinetic calculations
    bool use_smart_kinetic_solver = false;

    // Smart equilibrium method with clustering
    SmartKineticStrategy smart_kin_method = SmartKineticStrategy::Clustering;
    double smart_kinetic_tol = 1e-3;
    double smart_kinetic_abstol = 1e-4;
    double smart_kinetic_reltol = 1e-1;

    /// Create results file with parameters of the test
    auto makeResultsFile(const std::string& demo_tag) -> std::string
    {
        std::ostringstream eqtol_stream, kinreltol_stream, kinabstol_stream, kintol_stream;
        eqtol_stream << std::scientific << std::setprecision(1) << smart_equilibrium_reltol;
        kintol_stream << std::scientific << std::setprecision(1) << smart_kinetic_tol;
        kinreltol_stream << std::scientific << std::setprecision(1) << smart_kinetic_reltol;
        kinabstol_stream << std::scientific << std::setprecision(1) << smart_kinetic_abstol;

        std::string smart_test_tag = "-t0-" + std::to_string(int(t0)) +
                                     "-tfinal-" + std::to_string(int(tfinal)) +
                                     "-n-" + std::to_string(n) +
                                     "-eq-" + getSmartEquilibriumMethodTag() +
                                     "-kin-" + getSmartKineticMethodTag() +
                                     "-eqtol-" + eqtol_stream.str() +
                                     "-kintol-" + kintol_stream.str() +
                                     "-kinrel-" + kinreltol_stream.str() +
                                     "-kinabs-" + kinabstol_stream.str() +
                                     "-" + getActivityModelTag(activity_model) +
                                     (use_smart_kinetic_solver ? "-smart-kin" : "-conv-kin") +
                                     (use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq");      // name of the folder with results

        std::string class_test_tag = "-t0-" + std::to_string(int(t0)) +
                                     "-tfinal-" + std::to_string(int(tfinal)) +
                                     "-n-" + std::to_string(n) +
                                     "-" + getActivityModelTag(activity_model) +
                                     (use_smart_kinetic_solver ? "-smart-kin" : "-conv-kin") +
                                     (use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq");      // name of the folder with results
        std::string filename;
        if(use_smart_equilibrium_solver || use_smart_kinetic_solver)
            filename = "kineticpath-" + demo_tag + smart_test_tag + ".txt";
        else
            filename = "kineticpath-" + demo_tag + class_test_tag + ".txt";

        return filename;
    }
    auto getSmartEquilibriumMethodTag() -> std::string
    {
        switch(smart_eq_method)
        {
            case SmartEquilibriumStrategy::Clustering: return "clustering";
            case SmartEquilibriumStrategy::PriorityQueue: return "priority";
            case SmartEquilibriumStrategy::NearestNeighbour: return "nnsearch";
        }
        return "";
    }
    auto getSmartKineticMethodTag() -> std::string
    {
        switch(smart_kin_method)
        {
            case SmartKineticStrategy::Clustering: return "clustering";
            case SmartKineticStrategy::ClusteringExtended: return "clustering-extended";
            case SmartKineticStrategy::PriorityQueue: return "priority-major";
            case SmartKineticStrategy::PriorityQueuePrimary: return "priority-primary";
            case SmartKineticStrategy::NearestNeighbour: return "nnsearch";
        }
        return "";
    }

};

// Return string tag depending on the selected activity model
auto getActivityModelTag(enum ActivityModel activity_model) -> std::string
{
    switch(activity_model)
    {
        case ActivityModel::HKF: return "hkf-full";
        case ActivityModel::DebyeHuckel: return "dh-full";
        case ActivityModel::Pitzer: return "pitzer-full";
        case ActivityModel::HKFSelectedSpecies: return "hkf-selected-species";
        case ActivityModel::DebyeHuckelSelectedSpecies: return "dh-selected-species";
        case ActivityModel::PitzerSelectedSpecies: return "pitzer-selected-species";
    }
    return "";
}

// Return string tag depending on the selected method
auto getSmartMethodTag(enum SmartEquilibriumStrategy method) -> std::string
{
    switch(method)
    {
        case SmartEquilibriumStrategy::Clustering: return "eq-clustering";
        case SmartEquilibriumStrategy::PriorityQueue: return "eq-priority";
        case SmartEquilibriumStrategy::NearestNeighbour: return "eq-nnsearch";
    }
    return "";
}

// Return string tag depending on the selected method (smart kinetics strategy)
auto getSmartMethodTag(enum SmartKineticStrategy method) -> std::string
{
    switch(method)
    {
        case SmartKineticStrategy::Clustering: return "kin-clustering";
        case SmartKineticStrategy::PriorityQueue: return "kin-priority";
        case SmartKineticStrategy::NearestNeighbour: return "kin-nnsearch";
    }
    return "";
}

// Return string tag depending on the selected method fpr the hessian calculation
auto getGibbsHessianTag(enum GibbsHessian hessian) -> std::string
{
    switch(hessian)
    {
        case GibbsHessian::Exact: return "hessian-exact";
        case GibbsHessian::ExactDiagonal: return "hessian-exact-diagonal";
        case GibbsHessian::Approximation: return "hessian-approximation";
        case GibbsHessian::ApproximationDiagonal: return "hessian-approximation-diagonal";
    }
    return "";
}

/// Make directory for Windows and Linux
auto mkdir(std::string& folder) -> bool
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
