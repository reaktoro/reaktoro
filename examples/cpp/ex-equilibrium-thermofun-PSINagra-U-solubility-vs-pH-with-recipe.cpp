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
//   • Svetlana Kyas (27 September 2021)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <fstream>

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

    // Define chemical system by providing database, aqueous phase, and minerals
    ChemicalSystem system(db, solution);

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    // Create an equilibrium solver and its result
    EquilibriumSolver solver(specs);
    EquilibriumResult res;

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(25.0, "celsius");
    conditions.pressure(1.0, "bar");

    std::vector<double> pHs = {5, 6, 7, 8, 9, 10};
    std::vector<double> n(pHs);
    std::vector<double> nUO2(pHs);
    std::vector<double> nUO2OH(pHs);
    std::vector<double> nUO2OH2(pHs);
    std::vector<double> UO2CO3(pHs);        // UO2CO3@
    std::vector<double> UO23CO3OH3(pHs);    // (UO2)3CO3(OH)3+
    std::vector<double> UO22OH(pHs);        // (UO2)2(OH)+3
    std::vector<double> UO2CO33(pHs);       // UO2(CO3)3-4
    std::vector<double> UO2CO32(pHs);       // UO2(CO3)2-2
    std::vector<double> UO22OH2(pHs);       // (UO2)2(OH)2+2
    std::vector<double> UO23OH5(pHs);       // (UO2)3(OH)5+
    std::vector<double> UO24OH7(pHs);       // (UO2)4(OH)7+

    auto species_list = Strings{
        "UO2+2",
        "UO2OH+",
        "UO2(OH)2@",
        "UO2CO3@",
        "(UO2)3CO3(OH)3+",
        "(UO2)2(OH)+3",
        "UO2(CO3)3-4",
        "UO2(CO3)2-2",
        "(UO2)2(OH)2+2",
        "(UO2)3(OH)5+",
        "(UO2)4(OH)7+"
    };

    auto bU = 0.0;

    // ---------------------------------------------------------------------------------------------------------------//
    // OrdinaryPortlandCement-PoreWater
    // ---------------------------------------------------------------------------------------------------------------//

    // Result file for the portland cement problem
    std::ofstream result_file;
    result_file.precision(6);
    result_file.open("results-portlandcement.txt");

    // Define initial equilibrium state with the following recipe:
    // 1000 g H2O
    // H3PO4@ 1e-5 mol
    // CO2@ 0.1 g
    // UO2(SO4)@ 1e-5 mol.
    ChemicalState state(system);
    state.set("H2O@"     , 1e3,  "g");
    state.set("CO2@"     , 1e-1, "g");
    state.set("H3PO4@"   , 1e-5, "mol");
    state.set("UO2(SO4)@", 1e-5, "mol");

    // Aqueous properties of the chemical state
    AqueousProps aprops_pc(system);
    ChemicalProps props_pc(system);

    // Output the header of the table
    std::cout << "        pH       success         "
                 "%n(UO2+2)     %n(UO2OH+)   %n(UO2(OH)2@)    %n(UO2CO3@)   %n((UO2)3CO3(OH)3+)   "
                 "%n((UO2)2(OH)+3)   %n(UO2(CO3)3-4)   %n(UO2(CO3)2-2) %n((UO2)2(OH)2+2) %n((UO2)4(OH)7+)" << std::endl;

    for(Index i = 0; i < pHs.size(); i++)
    {
        conditions.pH(pHs[i]);

        // Equilibrate the initial state with given conditions and component amounts
        res = solver.solve(state, conditions);

        // Update aqueous properties
        aprops_pc.update(state);
        props_pc.update(state);

        //
        bU = props_pc.elementAmount("U");
        for(auto species : species_list)
        {
            result_file << state.speciesAmount(species) / bU * 100 << " ";
        }
        nUO2[i]       = state.speciesAmount("UO2+2");
        nUO2OH[i]     = state.speciesAmount("UO2OH+");
        nUO2OH2[i]    = state.speciesAmount("UO2(OH)2@");
        UO2CO3[i]     = state.speciesAmount("UO2CO3@");
        UO23CO3OH3[i] = state.speciesAmount("(UO2)3CO3(OH)3+");
        UO22OH[i]     = state.speciesAmount("(UO2)2(OH)+3");
        UO2CO33[i]    = state.speciesAmount("UO2(CO3)3-4");
        UO2CO32[i]    = state.speciesAmount("UO2(CO3)2-2");
        UO22OH2[i]    = state.speciesAmount("(UO2)2(OH)2+2");
        UO24OH7[i]    = state.speciesAmount("(UO2)4(OH)7+");

        // Output results for the current pH
        std::cout.precision(8);
        std::cout << std::scientific << " "
                  << pHs[i] << "    "
                  << res.succeeded() << "        "
                  << nUO2[i] / bU * 100 << " "
                  << nUO2OH[i] / bU * 100 << "  "
                  << nUO2OH2[i] / bU * 100 << " "
                  << UO2CO3[i] / bU * 100 << "        "
                  << UO23CO3OH3[i] / bU * 100 << "     "
                  << UO22OH[i] / bU * 100 << "    "
                  << UO2CO33[i] / bU * 100 << "    "
                  << UO2CO32[i] / bU * 100 << " "
                  << UO22OH2[i] / bU * 100 << " "
                  << UO24OH7[i] / bU * 100 << std::endl;

        // Output the chemical state to a console
        state.output("state-psinagra1207-portlandcement-with-recipe-pH-" + std::to_string(pHs[i]) + ".txt");
        result_file << "\n";
    }

    result_file.close();

    return 0;
}
