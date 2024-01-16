// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
//   • Svetlana Kyas (27 September August 2021)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Define Thermofun database
    ThermoFunDatabase db("mines16");

    // Define list of elements
    StringList selected_elements = "H O S Ca C Na Cl Fe Pb Zn";

    // Define aqueous phase
    AqueousPhase solution(speciate(selected_elements));
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    // Define minerals
    MineralPhases minerals("Arg Ang Gp Gn Sp Tro Smt Znc Ankerite Calcite "
                           "Siderite Graphite Fe(OH)3 Goethite Hematite "
                           "Lime Magnetite Anhydrite Pyrite Pyrrhotite Sulfur");

    // Define chemical system by providing database, aqueous phase, and minerals
    ChemicalSystem system(db, solution, minerals);

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.charge();
    specs.openTo("Cl-");

    // Create an equilibrium solver and its result
    EquilibriumSolver solver(specs);
    EquilibriumResult res;

    // Define initial state of the rock
    ChemicalState rockstate(system);
    //    H2SO4	    H2O	O2	        CaCO3	NaCl	    Fe	    Pb	        S	        Zn	        CH4
    //    1.00E-09	10	1.00E-09	1000	1.00E-09	0.002	1.00E-09	1.00E-09	1.00E-09	1.00E-09
    rockstate.setSpeciesMass("H2O@",      1e1, "g");
    rockstate.setSpeciesMass("Ca(CO3)@",  1e3, "g");
    rockstate.setSpeciesAmount("HSO4-",  1e-9, "mol");  // part of H2SO4
    rockstate.setSpeciesAmount("H+",     1e-9, "mol");  // part of H2SO4
    rockstate.setSpeciesAmount("O2@",    1e-9, "mol");
    rockstate.setSpeciesAmount("NaCl@",  1e-9, "mol");
    rockstate.setSpeciesAmount("Fe+3",   2e-3, "mol");
    rockstate.setSpeciesAmount("Pb+2",   1e-9, "mol");
    rockstate.setSpeciesAmount("Zn+2",   1e-9, "mol");
    rockstate.setSpeciesAmount("CH4@",   1e-9, "mol");
    rockstate.setSpeciesAmount("Sulfur", 1e-9, "mol");

    // Define temperature and pressure
    double T = 70.0; // in Celsius
    double P = 1.0; // in bar

    // Define conditions to be satisfied at chemical equilibrium of the rock
    EquilibriumConditions rockcond(specs);
    rockcond.temperature(T, "celsius");
    rockcond.pressure(P, "bar");
    rockcond.charge(0.0, "mol"); // to make sure the mixture is neutral

    // Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(rockstate, rockcond);
    std::cout << "res (rockstate) = " << res.succeeded() << std::endl;

    // Output the chemical state to a console
    rockstate.output("state-mines16-rock.txt");

    // Define initial state of the fluid
    ChemicalState fluidstate(system);
    //    H2SO4	    H2O	    O2	        CaCO3	    NaCl  Fe	    Pb	    S	    Zn	    CH4
    //    0.0556	1000	0.029499	1.00E-09	4	  0.002 	0.004	0.001	0.05	1.00E-09
    fluidstate.setSpeciesMass("H2O@",       1e3, "g");
    fluidstate.setSpeciesMass("Ca(CO3)@",  1e-9, "g");
    fluidstate.setSpeciesAmount("HSO4-", 0.0556, "mol"); // part of H2SO4
    fluidstate.setSpeciesAmount("H+",    0.0556, "mol"); // part of H2SO4
    fluidstate.setSpeciesAmount("O2@", 0.029499, "mol");
    fluidstate.setSpeciesAmount("NaCl@",    4.0, "mol");
    fluidstate.setSpeciesAmount("Fe+3",    2e-3, "mol");
    fluidstate.setSpeciesAmount("Pb+2",    4e-3, "mol");
    fluidstate.setSpeciesAmount("Zn+2",    5e-2, "mol");
    fluidstate.setSpeciesAmount("CH4@",    1e-9, "mol");
    fluidstate.setSpeciesAmount("Sulfur",  1e-3, "mol");

    // Define conditions to be satisfied at chemical equilibrium of the fluid
    EquilibriumConditions fluidcond(specs);
    fluidcond.temperature(T, "celsius");
    fluidcond.pressure(P, "bar");
    fluidcond.charge(0.0, "mol"); // to make sure the mixture is neutral

    // Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(fluidstate, fluidcond);
    std::cout << "res (fluidstate) = " << res.succeeded() << std::endl;

    // Output the chemical state to a console
    fluidstate.output("state-mines16-fluid.txt");

    return 0;
}
