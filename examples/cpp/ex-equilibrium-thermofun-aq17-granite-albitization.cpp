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

    // Define equilibrium solver
    EquilibriumOptions opts;

    EquilibriumSolver solver(system);

   // Define temperature and pressure
//    double T = 400.0; // in Celsius
//    double P = 1e3; // in bar

    double T = 25.0; // in Celsius
    double P = 1.0; // in bar

    std::cout << "T = " << T << std::endl;
    std::cout << "P = " << P << std::endl;

    // Initialize the amount of elements in the system
    Index E = system.elements().size();

    // Define initial equilibrium state for the granite calculations
    ChemicalState stategranite(system);
    stategranite.temperature(T, "celsius");
    stategranite.pressure(P, "bar");

    // -------------------------------------------------------------------------------------------- //
    // Pure granite
    // -------------------------------------------------------------------------------------------- //

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
    std::cout << "res (granite) = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    stategranite.output("state-aq17-granite.txt");

    // -------------------------------------------------------------------------------------------- //
    // Mix of granite and fluid
    // -------------------------------------------------------------------------------------------- //

    // Define initial equilibrium state for the granite-fluid mix calculations
    ChemicalState stategranitefluid(system);
    stategranitefluid.temperature(T, "celsius");
    stategranitefluid.pressure(P, "bar");

    // Define granite-fluid element amounts (mixed granit/fluid 0.2 mass ratio)
    // GEMS input:
    //Al e 0.84149656
    //Cl e 0.98929196
    //H h 104.59826
    //K e 0.2356788
    //Na e 1.3932998
    //O o 58.324214
    //Si e 2.2215454
    ArrayXr bgranitefluid(E + 1);
    // H, O, Na, Al, Si, Cl, K
    bgranitefluid << 104.59826, 58.324214, 1.3932998, 0.84149656, 2.2215454, 0.98929196, 0.2356788, 0.0;

    // Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(stategranitefluid, bgranitefluid);
    std::cout << "res (granite and fluid) = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    stategranitefluid.output("state-aq17-bgranitefluid.txt");

    // -------------------------------------------------------------------------------------------- //
    // Pure fluid
    // -------------------------------------------------------------------------------------------- //

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
    std::cout << "res (fluid) = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    statefluid.output("state-aq17-fluid.txt");

    return 0;
}