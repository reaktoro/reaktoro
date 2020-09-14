// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

int main()
{
    SupcrtDatabase db("supcrt98.xml");

    Phases phases(db,  AqueousSolution("H2O(l) H+ OH- O2(aq) H2(aq)"));

    ChemicalSystem system(db, phases);

    ChemicalState state(system);
    state.setTemperature(60.0, "celsius");
    state.setPressure(30.0, "bar");
    state.setSpeciesMass("H2O(l)", 1.0, "kg");

    EquilibriumOptions options;
    options.optima.output.active = true;

    EquilibriumSolver solver(system);
    solver.setOptions(options);

    solver.solve(state);

    std::cout << state.speciesAmounts() << std::endl;

    return 0;
}
