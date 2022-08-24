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
//   • Svetlana Kyas (27 September 2021)
//
// and since revised by:
//   • G. Dan Miron (28 January 2022)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Define Thermofun database
    ThermoFunDatabase db("aq17");

    // Define list of aqueous species
    StringList selected_species = "H2O@ H+ OH- Cl- HCl@ Na+ NaOH@ NaHSiO3@ NaCl@ NaAl(OH)4@ SiO2@ HSiO3- "
                                  "K+ KOH@ KCl@ KAlO2@ Al+3 AlOH+2 Al(OH)2+ Al(OH)3@ Al(OH)4-";

    // Define aqueous phase
    AqueousPhase solution(selected_species);

    // Set up a and b parameters for ionic species (NaCl, b = 0.064, a = 3.72)
    ActivityModelDebyeHuckelParams params;
    params.aiondefault = 3.72;
    params.biondefault = 0.064;
    params.bneutraldefault = 0.064;
    solution.setActivityModel(ActivityModelDebyeHuckel(params));

    // Define minerals
    MineralPhases minerals("Quartz Diaspore Gibbsite Andalusite Kyanite "
                           "Sillimanite Muscovite Paragonite Pyrophyllite "
                           "Kaolinite Albite Microcline");

    // Define chemical system by providing database, aqueous phase, and minerals
    ChemicalSystem system(db, solution, minerals);

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();

    // Create an equilibrium solver
    EquilibriumSolver solver(specs);

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(400.0, "celsius");
    conditions.pressure(1e3, "bar");

    // Define initial equilibrium state of 100 g of granite and 20 g of water
    ChemicalState state(system);

    // Define initial solution amount 1 NaCl m/Kg
    state.set("H2O@" , 20, "g");
    state.set("NaCl@", 0.02, "mol");

    // Define initial composition of granite (100g)
    state.set("Quartz"    , 35, "g"); // 35% of granite
    state.set("Microcline", 17, "g"); // 17% of granite
    state.set("Albite"    , 29, "g"); // 29% of granite
    state.set("Muscovite" , 19, "g"); // 19% of granite

    // Equilibrate the initial state with given conditions
    auto res = solver.solve(state, conditions);
    std::cout << "res (granite) = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    state.output("state-granite.txt");

    // Output the chemical properties to a text-file
    ChemicalProps props(state);
    props.output("props.txt");

    // Output the aqueous properties to a text-file
    AqueousProps aprops(state);
    aprops.output("aprops.txt");

    return 0;
}
