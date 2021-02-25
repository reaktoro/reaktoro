// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

/// **Note:**
/// This demo should be executed from the root directory of the binaries produced after compilation.
/// Example:
/// ~~~
/// cmake .
/// make -j
/// ./bin/demo-equilibrium-gems
/// ~~~
int main()
{
    // Use an exported project file from GEMS to initialize a Gems object,
    Gems gems("resources/gems/CalciteBC-dat.lst");

    // and then use it to construct the ChemicalSystem object.
    ChemicalSystem system = gems;

    // Create a ChemicalState object that contains the temperature, pressure,
    // and amounts of species stored in the exported GEMS file.
    ChemicalState state = gems.state(system);

    // Output the equilibrium state calculated by GEMS to a file.
    state.output("state-gems.txt");

    // Perturb the equilibrium state calculated by GEMS
    state.setSpeciesAmount("CO2@", 0.1);

    // and then equilibrate the modified chemical state using Reaktoro's methods.
    equilibrate(state);

    // Output the updated equilibrium state to a file.
    state.output("state-gems-updated.txt");
}
