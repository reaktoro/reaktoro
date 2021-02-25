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
    PhreeqcDatabase db("phreeqc.dat");

    AqueousPhase aqueousphase(speciate("H O C Na Cl"));
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gaseousphase("CO2(g)");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);

    ChemicalSystem system(phases);
    ChemicalState state(system);
    state.setTemperature(25.0, "celsius");
    state.setPressure(1.0, "bar");
    state.setSpeciesMass("H2O", 1.0, "kg");
    state.setSpeciesAmount("CO2(g)", 10.0, "mol");
    state.setSpeciesAmount("Na+", 4.0, "mol");
    state.setSpeciesAmount("Cl-", 4.0, "mol");

    EquilibriumOptions options;
    options.optima.output.active = true;

    EquilibriumSolver solver(system);
    solver.setOptions(options);

    solver.solve(state);

    const auto n = state.speciesAmounts();

    for(auto i = 0; i < n.size(); ++i)
    {
        std::cout << std::setw(20) << system.species(i).name();
        std::cout << std::setw(20) << n[i];
        std::cout << std::endl;
    }

    return 0;
}
