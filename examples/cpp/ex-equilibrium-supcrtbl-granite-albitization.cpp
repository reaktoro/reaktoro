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
//   • Svetlana Kyas (27 October 2021)
//
// and since revised by:
//   • Allan Leal (28 August 2023)
//     - Disabled the example by commenting out the code as it needs to be revised, with setup using Material class.
//   • Allan Leal (12 October 2022)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
//     // Initialize a thermodynamic database
//     SupcrtDatabase db("supcrtbl");

//     // Define list of aqueous species
//     StringList selected_species =
//         "H2O(aq) H+ OH- O2(aq) H2(aq) HCl(aq) Cl- SiO2(aq) HSiO3- "
//         "NaOH(aq) NaHSiO3(aq) NaCl(aq) NaAl(OH)4(aq) Na+ "
//         "KOH(aq) KCl(aq) KAlO2(aq) K+ "
//         "AlOH+2 Al+3 Al(OH)3(aq) Al(OH)4- Al(OH)2+";

//     // Define aqueous phase
//     AqueousPhase solution(selected_species);

//     // Set up a and b parameters for ionic species (NaCl, b = 0.064, a = 3.72)
//     ActivityModelDebyeHuckelParams params;
//     params.aiondefault = 3.72;
//     params.biondefault = 0.064;
//     params.bneutraldefault = 0.064;
//     solution.setActivityModel(ActivityModelDebyeHuckel(params));

//     // Define minerals
//     // Minerals that are not found: Cristobalite, Topaz-OH, Tridymite
//     MineralPhases minerals("Albite Andalusite Coesite Corundum Diaspore "
//         "Halite Kaolinite Kyanite Microcline Muscovite "
//         "Paragonite Pyrophyllite Quartz Sillimanite Stishovite "
//         "Sylvite");

//     // Define chemical system by providing database, aqueous phase, and minerals
//     ChemicalSystem system(db, solution, minerals);

//     // Set options for the equilibrium solver
//     EquilibriumOptions options;
//     // options.optima.output.active = true;

//     // Create an equilibrium solver
//     EquilibriumSolver solver(system);
//     solver.setOptions(options);

//     // The number of elements in the system
//     Index E = system.elements().size();

//     // The element amounts in the granite from GEMS
//     //   H  1.00E-16 mol
//     //   O  30.125894 mol
//     //   Na 2.0200391 mol
//     //   Al 4.2074828 mol
//     //   Si 11.107727 mol
//     //   Cl 1.00E-16 mol
//     //   K  1.178394 mol
//     ArrayXd bgranite(E + 1);
//     bgranite << 1.00e-16, 30.125894, 2.0200391, 4.2074828, 11.107727, 1.00e-16, 1.178394, 0.0; // H, O, Na, Al, Si, Cl, K, Z

//     // The element amounts in the fluid from GEMS
//     //   H  104.59826 mol
//     //   O  52.299035 mol
//     //   Na 0.98929196 mol
//     //   Al 1.00E-16 mol
//     //   Si 1.00E-16 mol
//     //   Cl 0.98929196 mol
//     //   K  1.00E-16 mol
//     ArrayXd bfluid(E + 1);
//     bfluid << 104.59826, 52.299035, 0.98929196, 1.00e-16, 1.00e-16, 0.98929196, 1.00e-16, 0.0; // H, O, Na, Al, Si, Cl, K, Z

//     // The element amounts in the entire system fluid and granite
//     ArrayXd b = bgranite + bfluid;

//     // The required conditions at equilbrium, including the amounts of components (elements and charge)
//     EquilibriumConditions conditions(system);
//     conditions.temperature(25.0, "celsius");
//     conditions.pressure(1.0, "bar");
//     conditions.setInitialComponentAmounts(b);

//     // The chemical state used to store the calculated equilibrium state below
//     ChemicalState state(system);

//     // Perform the equilibrium calculation with given element and charge amounts in `conditions`
//     options.optima.output.active = true;
//     solver.setOptions(options);

//     auto result = solver.solve(state, conditions);

//     // Check the did not fail
//     errorif(result.failed(), "The calculation did not succeed!");

//     // Output the computed chemical equilibrium state to a file
//     state.output("state-granite-fluid.txt");

//     return 0;
}
