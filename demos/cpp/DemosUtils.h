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

struct ReactiveTransportParams
{
    /// The options for the activity model of aqueous species
    enum struct AqueousActivityModel
    {
        HKF,

        DebyeHuckel,

        Pitzer,

        HKFSelectedSpecies,

        DebyeHuckelSelectedSpecies,

        PitzerSelectedSpecies,

    };

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

    AqueousActivityModel activity_model = AqueousActivityModel::HKF;

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
        std::cout << "activity model : " << getActivityModel() << std::endl;
        std::cout << "smart method   : " << getSmartMethodtag() << std::endl;

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
            case AqueousActivityModel::DebyeHuckel: return "dk-full";
            case AqueousActivityModel::Pitzer: return "pitzer-full";
            case AqueousActivityModel::HKFSelectedSpecies: return "hkf-selected-species";
            case AqueousActivityModel::DebyeHuckelSelectedSpecies: return "dk-selected-species";
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
                               "-" + getActivityModel() + "-reference";

        std::string smart_test_tag = "-" + getSmartMethodtag() +
                                     "-dt-" + dt_stream.str() +
                                     "-ncells-" + std::to_string(ncells) +
                                     "-nsteps-" + std::to_string(nsteps) +
                                     "-reltol-" + reltol_stream.str() +
                                     "-" + getActivityModel() + "-smart";

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

    // Discretisation params
    double t0;      // starting time of simulations
    double tfinal;  // final time of simulation

    /// Create results file with parameters of the test
    /// Create results file with parameters of the test
    auto makeResultsFolder() -> std::string
    {
        struct stat status = {0};               // structure to get the file status

        std::ostringstream tol_stream, t0_stream, tfinal_stream;
        t0_stream << t0;
        tfinal_stream << tfinal;

        std::string test_tag = "-t0-" + t0_stream.str() +
                               "-tfinal-" + tfinal_stream.str();      // name of the folder with results
        std::string folder = "../kinetics-perturbed-co2" + test_tag;
        if (stat(folder.c_str(), &status) == -1) mkdir(folder);

        return folder;
    }

};
