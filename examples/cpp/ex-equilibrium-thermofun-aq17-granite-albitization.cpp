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
    StringList selected_species =
        "H2O@ H+ OH- O2@ H2@ HCl@ Cl- SiO2@ HSiO3- "
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

    // Create an equilibrium solver
    EquilibriumOptions opts;

    // Create equilibrium solver
    EquilibriumSolver solver(system);

    // Create the result of the equilibrium calculations
    EquilibriumResult res;

    // Initialize the amount of elements in the system
    Index E = system.elements().size();

    // Define initial equilibrium state for the equilibrium calculations
    ChemicalState initialstate(system);
    initialstate.temperature(400.0, "celsius");
    initialstate.pressure(1e3, "bar");
    initialstate.set("H2O@", 1.0, "kg");

    // -------------------------------------------------------------------------------------------- //
    // Mix of granite and fluid
    // -------------------------------------------------------------------------------------------- //

    // Define initial equilibrium state for the granite-fluid mix calculations
    ChemicalState stategranitefluid = initialstate;

    // Define granite-fluid element amounts (mixed granit/fluid 0.2 mass ratio)
    // GEMS input:
    // Al e 0.84149656
    // Cl e 0.98929196
    // H h 104.59826
    // K e 0.2356788
    // Na e 1.3932998
    // O o 58.324214
    // Si e 2.2215454
    // Z 0.0
    ArrayXd bgranitefluid(E + 1);
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
    ChemicalState statefluid = initialstate;

    // Define granite element amounts
    // GEMS input:
    //    Al    e   	1.00E-09
    //    Cl    e   	0.98929196
    //    H     h   	104.59826
    //    K     e   	1.00E-09
    //    Na    e   	0.98929196
    //    O     o   	52.299035
    //    Si    e   	1.00E-09
    //    Z             0.0
    ArrayXd bfluid(E + 1);
    // H, O, Na, Al, Si, Cl, K
    bfluid << 104.59826, 52.299035, 0.98929196, 1.00e-09, 1.00e-09, 0.98929196, 1.00e-09, 0.0;

    // Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(statefluid, bfluid);
    std::cout << "res (fluid) = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    statefluid.output("state-aq17-fluid.txt");

    return 0;
}
