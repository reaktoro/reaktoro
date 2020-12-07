// Reaktoro is a unified framework for modeling chemically reactive systems.
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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

// This string defines a PHREEQC script problem.
// This problem was taken from the official PHREEQC example named ex1.
const std::string ex3 = R"(
TITLE Example 3, part A.--Calcite equilibrium at log Pco2 = -2.0 and 25C.
SOLUTION 1  Pure water
        pH      7.0
        temp    25.0
EQUILIBRIUM_PHASES
        CO2(g)          -2.0
        Calcite         0.0
SAVE solution 1
END
TITLE Example 3, part B.--Definition of seawater.
SOLUTION 2  Seawater
        units   ppm
        pH      8.22
        pe      8.451
        density 1.023
        temp    25.0
        Ca              412.3
        Mg              1291.8
        Na              10768.0
        K               399.1
        Si              4.28
        Cl              19353.0
        Alkalinity      141.682 as HCO3
        S(6)            2712.0
END
TITLE Example 3, part C.--Mix 70% groundwater, 30% seawater.
MIX 1
        1      0.7
        2      0.3
SAVE solution   3
END
TITLE Example 3, part E.--Equilibrate mixture with calcite only.
EQUILIBRIUM_PHASES 2
        Calcite         0.0
USE solution 3
END
)";

int main()
{
    ChemicalEditor editor;
    editor.addAqueousPhaseWithElements("H O C Ca Cl Na K Mg S Si");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Dolomite");

    ChemicalSystem system(editor);

    // Amount of water
    double water_kg = 1.00;
    double T = 25.0;
    double P = 1.0;

    // Define the state with pure water
    EquilibriumInverseProblem problem_pw(system);
    problem_pw.setTemperature(T, "celsius");
    problem_pw.setPressure(P, "atm");
    problem_pw.add("H2O", water_kg, "kg");
    problem_pw.pH(7.0, "HCl", "NaOH");

    ChemicalState state_pw = equilibrate(problem_pw);
    std::cout << "Pure water: \n" << state_pw << std::endl;

    // Define the state with seawater
    EquilibriumInverseProblem problem_sw(system);
    problem_sw.setTemperature(T, "celsius");
    problem_sw.setPressure(P, "atm");
    problem_sw.add("H2O", water_kg, "kg");
    problem_sw.add("SO4--", 2712.0 * water_kg, "mg");
    problem_sw.add("Ca++", 412.3 * water_kg, "mg");
    problem_sw.add("Na+", 10768.0 * water_kg, "mg");
    problem_sw.add("K+", 399.1 * water_kg, "mg");
    problem_sw.add("Mg++", 1290 * water_kg, "mg");
    problem_sw.add("Cl-", 19353.0 * water_kg, "mg");
    problem_sw.add("HCO3-", 141.682 * water_kg, "mg");
    problem_sw.add("Si", 4.28 * water_kg, "mg");
    problem_sw.pH(8.22, "HCl", "NaOH");

    ChemicalState state_sw = equilibrate(problem_sw);
    std::cout << "Sea water: \n" << state_sw << std::endl;

    // Conduct the mixing pure water and seawater
    ChemicalState state_mixed = 0.3 * state_sw + 0.7 * state_pw;
    state_mixed.setSpeciesAmount("Calcite", 10.0, "mol");

    // Define solver for equilibration new mixed state
    EquilibriumSolver solver(system);
    EquilibriumResult res = solver.solve(state_mixed);
    std::cout << "Mixture(SW 30% and PW 70%) equilibrated with 10 mol of Calcite: \n" << state_mixed << std::endl;

    // Initialize a Phreeqc instance with the official phreeqc.dat database file
    Phreeqc phreeqc("databases/phreeqc/phreeqc.dat");

    // Execute a PHREEQC script defining a geochemical problem.
    // Here this script is actually embedded into a string named `ex1`.
    // However, `ex3` could also be a string containing the path to a script file.
    // Method execute will automatically identify when the contents are embedded in the string and
    // when the string is actually a path to a script file.
    phreeqc.execute(ex3);

    // Initialize a ChemicalSystem instance using the current state of the Phreeqc instance.
    // This will allow the use of both PHREEQC thermodynamic data and PHREEQC activity models
    // in the subsequent equilibrium calculations using Reaktoro's algorithms.
    ChemicalSystem system_phreeqc(phreeqc);

    // Initialize an ChemicalState instance using the current state of the Phreeqc instance.
    ChemicalState state_phreeqc = phreeqc.state(system_phreeqc);

    // Print such equilibrium state using Reaktoro format.
    std::cout << "PHREEQC mix with 70% PW and 30% SW equilibrated with Calcite: \n " << state_phreeqc << std::endl;
}
