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
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Define Thermofun database
    ThermoFunDatabase db("aq17");

    // Define list of aqueous species
    StringList selected_species = "H2O@ H+ OH- Cl- HCl@ Na+ NaOH@ NaHSiO3@ NaCl@ NaAl(OH)4@ "
                                  "K+ KOH@ KCl@ KAlO2@ Al+3 AlOH+2 Al(OH)2+ Al(OH)3@ Al(OH)4-";

    // Define aqueous phase
    AqueousPhase solution(selected_species);
    solution.setActivityModel(ActivityModelHKF());

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

    // Define equilibrium solver
    EquilibriumSolver solver(specs);

    // Define temperature and pressure
    double T = 300.0; // in Celsius
    double P = Reaktoro::waterSaturatedPressureWagnerPruss(T + 273.15).val() * 1e-5; // is Psat of water at the T = 300

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(T, "celsius");
    conditions.pressure(P, "bar");

    // Define initial solution amount
    conditions.startWith("H2O@", 55.51, "mol");
    conditions.startWith("NaCl@", 0.27, "mol");
    conditions.startWith("KCl@", 0.03, "mol");

    // Define initial granite composition
    conditions.startWith("Quartz"    , 168.126, "mol"); //
    conditions.startWith("Microcline", 17.8099, "mol"); // K-Feldspar, K(AlSi3)O8 + 8 * H2O = K+ + Al(OH)4- + 6 * H2O + 3 * SiO2(aq)
    conditions.startWith("Albite"    , 19.937 , "mol"); // Na(AlSi3)O8 + 8 * H2O = Na+ + Al(OH)4- + 6 * H2O + 3 * SiO2(aq)
    conditions.startWith("Muscovite" , 2.15255, "mol"); // KAl2(AlSi3)O10(OH)2 = 1K+ + 3 * Al+++ - 10 * H+ + 6 * H2O + 3 * SiO2(aq)

    // Define initial equilibrium state
    ChemicalState state(system);

    // Equilibrate the initial state with given conditions
    solver.solve(state, conditions);

    // Output the chemical state to a console
    state.output("state.txt");

    return 0;
}
