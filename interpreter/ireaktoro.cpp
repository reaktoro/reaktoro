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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

/// Prints a message detailing how to use the interpreter executable.
void printUsage(int argc, char** argv);

/// The entry point of the Reaktoro interpreter
int main(int argc, char** argv)
{
    // Check the command-line arguments were used correctly
    if(argc != 2) {
        printUsage(argc, argv);
        return 0;
    }

    // Create the interpreter to execute the input file
    Interpreter interpreter;

    // Execute the input file provided via command-line
    interpreter.executeJsonFile(argv[1]);

    // Output all calculated states
    for(const auto& pair : interpreter.states())
        pair.second.output(pair.first + ".txt");
}

void printUsage(int argc, char** argv)
{
    std::cout << "Usage: " << argv[0] << " `your-input-file.json`" << std::endl;
}
