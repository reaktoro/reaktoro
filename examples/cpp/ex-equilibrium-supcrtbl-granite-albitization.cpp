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
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Define Reaktoro database
    SupcrtDatabase db("supcrtbl");

//    for(auto s : db.species()) {
//        std::cout << s.name() << std::endl;
//    }
//    getchar();
    // Define list of aqueous species
    StringList selected_species = "H2O(aq) H+ OH- O2(aq) H2(aq) HCl(aq) Cl- SiO2(aq) HSiO3- "
                                  "NaOH(aq) NaHSiO3(aq) NaCl(aq) NaAl(OH)4(aq) Na+ "
                                  "KOH(aq) KCl(aq) KAlO2(aq) K+ "
                                  "AlOH+2 Al+3 Al(OH)3(aq) Al(OH)4- Al(OH)2+";

    // Define aqueous phase
    AqueousPhase solution(selected_species);

    // Set up a and b parameters for ionic species (NaCl, b = 0.064, a = 3.72)
    ActivityModelDebyeHuckelParams params;
    params.aiondefault = 3.72;
    params.biondefault = 0.064;
    params.bneutraldefault = 0.064;
    solution.setActivityModel(ActivityModelDebyeHuckel(params));

    // Define minerals
    // Minerals that are not found: Cristobalite, Topaz-OH, Tridymite
    MineralPhases minerals("Albite Andalusite Coesite Corundum Diaspore "
                           "Halite Kaolinite Kyanite Microcline Muscovite "
                           "Paragonite Pyrophyllite Quartz Sillimanite Stishovite "
                           "Sylvite");

    // Define chemical system by providing database, aqueous phase, and minerals
    ChemicalSystem system(db, solution, minerals);

    // Specify conditions to be satisfied at chemical equilibrium
//    EquilibriumSpecs specs(system);
//    specs.temperature();
//    specs.pressure();

    // Define equilibrium solver
    EquilibriumOptions opts;

    EquilibriumSolver solver(system);

//    // Define temperature and pressure
//    double T = 400.0; // in Celsius
//    double P = 1e3; // in bar

    double T = 25.0; // in Celsius
    double P = 1.0; // in bar

    // Define initial equilibrium state for the granite calculations
    ChemicalState stategranite(system);
    stategranite.temperature(T, "celsius");
    stategranite.pressure(P, "bar");

    // Initialize the amount of elements in the system
    Index E = system.elements().size();

//    // Define conditions to be satisfied at chemical equilibrium
//    EquilibriumConditions conditions(specs);
//    conditions.temperature(T, "celsius");
//    conditions.pressure(P, "bar");

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
    opts.optima.output.active = false;
    solver.setOptions(opts);
    auto res = solver.solve(stategranite, bgranite);
    std::cout << "res = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    stategranite.output("state-granite.txt");

    // Define initial equilibrium state for the fluid calculations
    ChemicalState statefluid(system);
    statefluid.temperature(T, "celsius");
    statefluid.pressure(P, "bar");

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
    opts.optima.output.active = false;
    solver.setOptions(opts);
    res = solver.solve(statefluid, bfluid);
    std::cout << "res = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    statefluid.output("state-fluid.txt");

    return 0;
}
