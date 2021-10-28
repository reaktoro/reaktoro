// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Allan Leal (19 July 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution(speciate("H O C Na Cl Ca Mg Si"), exclude("organic"));
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gases("CO2(g) H2O(g)");
    gases.setActivityModel(ActivityModelPengRobinson());

    MineralPhases minerals;

    ChemicalSystem system(db, solution, gases, minerals);

    ChemicalState state(system);
    state.temperature(60.0, "celsius");
    state.pressure(100.0, "bar");
    state.set("H2O(aq)"  , 1.0, "kg");
    state.set("CO2(g)"   , 1.0, "mol");
    state.set("Halite"   , 1.0, "mol");
    state.set("Calcite"  , 1.0, "mol");
    state.set("Magnesite", 1.0, "mol");
    state.set("Quartz"   , 1.0, "mol");

    EquilibriumSolver solver(system);
    solver.solve(state);

    ChemicalProps props(state);
    props.output("props.txt");

    AqueousProps aprops(state);
    aprops.output("aprops.txt");

    std::cout << "Success! Check outputted files `props.txt` and `aprops.txt`." << std::endl;

    return 0;
}
