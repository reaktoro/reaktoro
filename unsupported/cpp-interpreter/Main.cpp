// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

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
