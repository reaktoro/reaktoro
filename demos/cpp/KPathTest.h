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

struct Params{

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
