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

int main()
{
    ThermoFunDatabase db("aq17");

    AqueousPhase aqueousphase(speciate("H O C Na Cl"));
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gaseousphase("CO2 H2O");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);

    ChemicalSystem system(phases);

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    EquilibriumSolver solver(specs);

    EquilibriumConditions conditions(specs);
    conditions.temperature(60.0, "celsius");
    conditions.pressure(100.0, "bar");
    conditions.pH(4.0);
    conditions.startWith("H2O@", 1.0, "kg");
    conditions.startWith("Na+",  1.0, "mol");
    conditions.startWith("Cl-",  1.0, "mol");
    conditions.startWith("CO2", 10.0, "mol");

    ChemicalState state(system);

    solver.solve(state, conditions);

    ChemicalProps props(state);
    AqueousProps aprops(state);

    std::cout << state << std::endl;
    std::cout << props << std::endl;
    std::cout << aprops << std::endl;

    return 0;
}
