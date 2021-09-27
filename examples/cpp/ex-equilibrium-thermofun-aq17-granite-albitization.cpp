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
    StringList selected_species = "H2O@ H+ OH- O2@ H2@ HCl@ Cl- SiO2@ HSiO3- "
                                  "NaOH@ NaHSiO3@ NaCl@ NaAl(OH)4@ Na+ "
                                  "KOH@ KCl@ KAlO2@ K+ "
                                  "AlOH+2 Al+3 Al(OH)3@ Al(OH)4- Al(OH)2+";

    // Define aqueous phase
    AqueousPhase solution(selected_species);

    // Set up a and b parameters for ionic species (NaCl, b = 0.064, a = 3.72)
    ActivityModelDebyeHuckelParams params;
    params.aiondefault = 3.72;
    params.biondefault = 0.064;
    params.bneutraldefault = 0.064;
    solution.setActivityModel(ActivityModelDebyeHuckel(params));

    // Define minerals
    MineralPhases minerals("Albite Andalusite Coesite Corundum Cristobalite Diaspore "
                           "Halite Kaolinite Kyanite Microcline Muscovite "
                           "Paragonite Pyrophyllite Quartz Sillimanite Stishovite "
                           "Sylvite Topaz-OH Tridymite");

    // Define chemical system by providing database, aqueous phase, and minerals
    ChemicalSystem system(db, solution, minerals);

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();

    // Define equilibrium solver
    EquilibriumSolver solver(specs);

    // Define initial equilibrium state
    ChemicalState state(system);

    // Define temperature and pressure
    double T = 400.0; // in Celsius
    double P = 1e3; // in bar

    // Initialize the amount of elements in the system
    Index E = system.elements().size();

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(T, "celsius");
    conditions.pressure(P, "bar");

    // Define granite element amounts
    // GEMS input:
    //    Al    e   	4.2074828
    //    Cl    e   	1.00E-09
    //    H     h   	1.00E-09
    //    K     e   	1.178394
    //    Na    e   	2.0200391
    //    O     o   	30.125894
    //    Si    e   	11.107727
    ArrayXr bgranite(E + 1);
    // H, O, Na, Al, Si, Cl, K
    bgranite << 1.00e-09, 30.125894, 2.0200391, 4.2074828, 11.107727, 1.00e-09, 1.178394, 0.0;

    // Equilibrate the initial state with given conditions and component amounts
    solver.solve(state, conditions, bgranite);

    // Output the chemical state to a console
    state.output("state-aq17-granite.txt");

    // Define granite element amounts
    // GEMS input:
    //    Al    e   	1.00E-09
    //    Cl    e   	0.98929196
    //    H     h   	104.59826
    //    K     e   	1.00E-09
    //    Na    e   	0.98929196
    //    O     o   	52.299035
    //    Si    e   	1.00E-09
    ArrayXr bfluid(E + 1);
    // H, O, Na, Al, Si, Cl, K
    bfluid << 104.59826, 52.299035, 0.98929196, 1.00e-09, 1.00e-09, 0.98929196, 1.00e-09, 0.0;

    // Equilibrate the initial state with given conditions and component amounts
    solver.solve(state, conditions, bfluid);

    // Output the chemical state to a console
    state.output("state-aq17-fluid.txt");

    return 0;
}