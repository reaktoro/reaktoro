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
    ThermoFunDatabase db("psinagra-12-07");

    // Define list of elements
    StringList selected_elements = "Al C Ca Cl Fe H K Mg N Na O P S Si U";

    // Define aqueous phase
    AqueousPhase solution(speciate(selected_elements));
    solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond("CO2")
    ));

    // Define minerals
    MineralPhases minerals(speciate(selected_elements));

    // Define chemical system by providing database, aqueous phase, and minerals
    ChemicalSystem system(db, solution, minerals);

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    // Define equilibrium solver and its result
    EquilibriumSolver solver(specs);
    EquilibriumResult res;

    // Define temperature and pressure
    double T = 25.0; // in Celsius
    double P = 1.0; // in bar

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(T, "celsius");
    conditions.pressure(P, "bar");
    conditions.pH(6);

    // Define initial equilibrium state
    ChemicalState solutionstate_pc(system);
    solutionstate_pc.setSpeciesMass("H2O@",   1e3,  "g");
    solutionstate_pc.setSpeciesMass("CO2@",   1e-1, "g");
    solutionstate_pc.setSpeciesMass("H3PO4@", 1e-5, "g");
    solutionstate_pc.setSpeciesMass("SO2@",   1e-1, "g");

    // Define element amounts for OrdinaryPortlandCement-PoreWater
    // GEMS input (OrdinaryPortlandCement-PoreWater 1 kg):
    //    Al    e   	8.06E-06
    //    C     e   	4.01E-05
    //    Ca    e   	2.53E-03
    //    Cl    e   	7.98E-08
    //    Fe    e   	3.28E-08
    //    H     h   	110.1607
    //    K     e   	0.13369924
    //    Mg    e   	4.60E-09
    //    N     e   	7.98E-08
    //    Na    e   	0.038244662
    //    O     o   	55.17178
    //    P     e   	7.98E-08
    //    S     e   	7.61E-04
    //    Si    e   	1.78E-05
    ArrayXr portlandcementb(system.elements().size() + 1);
    // H C N O Na Mg Al Si P S Cl K Ca Fe U Z
    portlandcementb << 110.1607, 4.01E-05, 7.98E-08, 55.17178, 0.038244662, 4.60E-09, 8.06E-06,
                        1.78E-05, 7.98E-08, 7.61E-04, 7.98E-08, 0.13369924, 2.53E-03, 3.28E-08,
                        1e-5, 0.0;

    // Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(solutionstate_pc, conditions, portlandcementb);
    std::cout << "res (OrdinaryPortlandCement-PoreWater) = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    solutionstate_pc.output("state-psinagra1207-portlandcement-pH6.txt");

    // Output the properties of the chemical state
    AqueousProps aprops_pc(solutionstate_pc);
    aprops_pc.output("aprops-psinagra1207-portlandcement-pH6.txt");

    // Define initial equilibrium state
    ChemicalState solutionstate_bc(system);
    solutionstate_bc.setSpeciesMass("H2O@",   1e3,  "g");
    solutionstate_bc.setSpeciesMass("CO2@",   1e-1, "g");
    solutionstate_bc.setSpeciesMass("H3PO4@", 1e-5, "g");
    solutionstate_bc.setSpeciesMass("SO2@",   1e-1, "g");

    // Define element amounts for Boda-Caly-pore-water 1 kg
    // GEMS input (Boda-Caly-pore-water 1 kg 1 kg):
    //    C     e   	0.00061
    //    Ca    e   	0.00311
    //    Cl    e   	0.023787
    //    H     h   	110.75091
    //    K     e   	0.00018
    //    Mg    e   	0.0024
    //    Na    e   	0.017
    //    O     o   	55.384583
    //    S     e   	0.0019
    ArrayXr bodycalyb(system.elements().size() + 1);
    // H C N O Na Mg Al Si P S Cl K Ca Fe U Z
    bodycalyb << 110.75091, 0.00061, 0.0, 55.384583, 0.017, 0.0024,
                0.0, 0.0, 0.0, 0.0019, 0.023787, 0.00018, 0.00311,
                0.0, 1e-5, 0.0;

    // Equilibrate the initial state with given conditions and component amounts
    res = solver.solve(solutionstate_bc, conditions, bodycalyb);
    std::cout << "res (Boda-Caly-pore-water) = " << res.optima.succeeded << std::endl;

    // Output the chemical state to a console
    solutionstate_bc.output("state-psinagra1207-bodycaly-pH6.txt");

    // Output the properties of the chemical state
    AqueousProps aprops_bc(solutionstate_bc);
    aprops_bc.output("aprops-psinagra1207-bodycaly-pH6.txt");

    return 0;
}
