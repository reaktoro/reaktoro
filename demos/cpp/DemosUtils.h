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


/// Make directory for Windows and Linux
auto mkdir(std::string& folder) -> bool
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
auto getActivityModelTag(enum ActivityModel activity_model) -> std::string;
auto mkdir(const std::string& folder) -> bool;


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
    double smart_equilibrium_reltol = 0;
    double amount_fraction_cutoff = 0;
    double mole_fraction_cutoff = 0;

    ActivityModel activity_model = ActivityModel::HKF;

    SmartEquilibriumStrategy method = SmartEquilibriumStrategy::Clustering;

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
        std::cout << "eqreltol       : " << smart_equilibrium_reltol << std::endl;
        std::cout << "activity model : " << getActivityModelTag(activity_model) << std::endl;
        std::cout << "smart method   : " << getSmartMethodTag(method) << std::endl;

    }

    auto getSmartMethodtag() -> std::string
    {
        switch(method)
        {
            case SmartEquilibriumStrategy::Clustering: return "eq-clustering";
            case SmartEquilibriumStrategy::PriorityQueue: return "eq-priority";
            case SmartEquilibriumStrategy::NearestNeighbour: return "eq-nnsearch";
        }
    }

    auto getActivityModel() -> std::string
    {
        switch(activity_model)
        {
            case AqueousActivityModel::HKF: return "hkf-full";
            case AqueousActivityModel::DebyeHuckel: return "dh-full";
            case AqueousActivityModel::Pitzer: return "pitzer-full";
            case AqueousActivityModel::HKFSelectedSpecies: return "hkf-selected-species";
            case AqueousActivityModel::DebyeHuckelSelectedSpecies: return "dh-selected-species";
            case AqueousActivityModel::PitzerSelectedSpecies: return "pitzer-selected-species";
        }
    }

    /// Create results file with parameters of the test
    auto makeResultsFolder(const std::string& demo_tag) -> std::string
    {
        struct stat status = {0};               // structure to get the file status

        std::ostringstream reltol_stream, dt_stream;
        dt_stream << dt;
        reltol_stream << std::scientific << std::setprecision(1) << smart_equilibrium_reltol;

        std::string test_tag = "-dt-" + dt_stream.str() +
                               "-ncells-" + std::to_string(ncells) +
                               "-nsteps-" + std::to_string(nsteps) +
                               "-" + getActivityModelTag(activity_model) + "-reference";

        std::string smart_test_tag = "-" + getActivityModelTag(activity_model) +
                                     "-dt-" + dt_stream.str() +
                                     "-ncells-" + std::to_string(ncells) +
                                     "-nsteps-" + std::to_string(nsteps) +
                                     "-reltol-" + reltol_stream.str() +
                                     "-" + getActivityModelTag(activity_model) + "-smart";

        std::string folder = "results-" + demo_tag;
        folder = (use_smart_equilibrium_solver) ?
                 folder + smart_test_tag :
                 folder + test_tag;

        if (stat(folder.c_str(), &status) == -1) mkdir(folder);

        std::cout << "\nsolver                         : " << (use_smart_equilibrium_solver ? "smart" : "conventional") << std::endl;

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

struct KineticPathParams{

    // Discretization params
    double t0 = 0;      // starting time of simulations
    double tfinal = 0;  // final time of simulation

    double dt = 0.0;    // discretization step in time
    int n = 0;          // number of steps

    // Thermodynamic parameters
    double T = 0; // the temperature (in units of degC)
    double P = 0; // the pressure (in units of bar)

    ActivityModel activity_model = ActivityModel::HKF;

    // ODML parameters
    bool use_smart_equilibrium_solver = false;
    double smart_equilibrium_reltol = 1e-3;

    SmartEquilibriumStrategy method = SmartEquilibriumStrategy::Clustering;

    /// Create results file with parameters of the test
    auto makeResultsFile(const std::string& demo_tag) -> std::string
    {
        std::ostringstream eqtol_stream, kinreltol_stream, kinabstol_stream, kintol_stream;
        eqtol_stream << std::scientific << std::setprecision(1) << smart_equilibrium_reltol;
        //kinreltol_stream << std::scientific << std::setprecision(1) << smart_kinetic_options.reltol;
        //kinabstol_stream << std::scientific << std::setprecision(1) << smart_kinetic_options.abstol;
        //kintol_stream << std::scientific << std::setprecision(1) << smart_kinetic_options.tol;

        std::string smart_test_tag = "-t0-" + std::to_string(t0) +
                                     "-tfinal-" + std::to_string(tfinal) +
                                     "-n-" + std::to_string(n) +
                                     "-" + getSmartMethodTag() +
                                     "-eqtol-" + eqtol_stream.str() +
                                     //"-kintol-" + kintol_stream.str() +
                                     //"-kinrel-" + kinreltol_stream.str() +
                                     //"-kinabs-" + kinabstol_stream.str() +
                                     "-" + getActivityModelTag(activity_model) +
                                     //(use_smart_kinetic_solver ? "-smart-kin" : "-conv-kin") +
                                     (use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq");      // name of the folder with results

        std::string class_test_tag = "-t0-" + std::to_string(t0) +
                                     "-tfinal-" + std::to_string(tfinal) +
                                     "-n-" + std::to_string(n) +
                                     "-" + getActivityModelTag(activity_model) +
                                     //(use_smart_kinetic_solver ? "-smart-kin" : "-conv-kin") +
                                     (use_smart_equilibrium_solver ? "-smart-eq"  : "-conv-eq");      // name of the folder with results
        std::string filename;
        if(use_smart_equilibrium_solver)
            filename = "kineticpath-" + demo_tag + smart_test_tag + ".txt";
        else
            filename = "kineticpath-" + demo_tag + class_test_tag + ".txt";

        return filename;
    }
    auto getSmartMethodTag() -> std::string
    {
        switch(method)
        {
            case SmartEquilibriumStrategy::Clustering: return "clustering";
            case SmartEquilibriumStrategy::PriorityQueue: return "priority";
            case SmartEquilibriumStrategy::NearestNeighbour: return "nnsearch";
        }
    }

};

// Return string tag depending on the selected activity model
auto getActivityModelTag(enum ActivityModel activity_model) -> std::string
{
    switch(activity_model)
    {
        case ActivityModel::HKF: return "hkf-full";
        case ActivityModel::DebyeHuckel: return "dk-full";
        case ActivityModel::Pitzer: return "pitzer-full";
        case ActivityModel::HKFSelectedSpecies: return "hkf-selected-species";
        case ActivityModel::DebyeHuckelSelectedSpecies: return "dk-selected-species";
        case ActivityModel::PitzerSelectedSpecies: return "pitzer-selected-species";
    }
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
