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
#include <unsupported/cpp-interpreter/Interpreter.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Reaktoro includes
using namespace Reaktoro;

int main(int argc, char **argv)
{
    // Collect the command-line arguments into a list of strings
    std::vector<std::string> args;
    for(int i = 0; i < argc; ++i)
        args.push_back(argv[i]);

    if(argc < 2)
    {
        std::cout << "Usage: reaktoro `scriptfile`" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::ifstream inputscript(filename);

    Interpreter interpreter;
    interpreter.execute(inputscript);

    return 0;
}
