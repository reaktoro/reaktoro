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
enum RTActivityModel
{
    HKF,

    DebyeHuckel,

    Pitzer,

};

// Forward declaration
auto getRTActivityModelTag(enum RTActivityModel activity_model) -> std::string;
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

    RTActivityModel activity_model = RTActivityModel::HKF;

    GibbsHessian hessian = GibbsHessian::Exact;

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
        std::cout << "activity model  : " << getRTActivityModelTag(activity_model) << std::endl;
        std::cout << "hessian         : " << getGibbsHessianTag(hessian) << std::endl;

    }

    /// Output test parameters to the console for the reactive transport with equilibrium only
    auto outputEquilibriumMethod() -> void
    {
        std::cout << "smart eq used   : " << use_smart_equilibrium_solver << std::endl;
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
                               "-" + getRTActivityModelTag(activity_model) + "-conventional";

        std::string folder = "rt-" + demo_tag + test_tag;

        if (stat(folder.c_str(), &status) == -1) mkdir(folder);

        //std::cout << "\nsolver                         : " << (use_smart_solver ? "smart" : "conventional") << std::endl;

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

    // Rate of the smart equilibrium estimation w.r.t to the total chemical equilibrium calculation
    double smart_equilibrium_acceptance_rate = 0.0;
};
// Return string tag depending on the selected activity model
auto getRTActivityModelTag(enum RTActivityModel activity_model) -> std::string
{
    switch(activity_model)
    {
        case RTActivityModel::HKF: return "hkf";
        case RTActivityModel::DebyeHuckel: return "dh";
        case RTActivityModel::Pitzer: return "pitzer";
    }
    return "";
}

// Return string tag depending on the selected method fpr the hessian calculation
auto getGibbsHessianTag(enum GibbsHessian hessian) -> std::string
{
    switch(hessian)
    {
        case GibbsHessian::Exact: return "hessian-exact";
        case GibbsHessian::PartiallyExact: return "hessian-partially-exact";
        case GibbsHessian::Approx: return "hessian-approx";
        case GibbsHessian::ApproxDiagonal: return "hessian-approx-diagonal";
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