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

//    // Define minerals
//    MineralPhases minerals(speciate(selected_elements));

    // Define chemical system by providing database, aqueous phase, and minerals
    //ChemicalSystem system(db, solution, minerals);
    ChemicalSystem system(db, solution);

//    //auto species_with_u = system.species().withElements("U");
//    auto species_with_u = system.species()
//            .withSomeOfElements("U")
//            .withAggregateState(AggregateState::Aqueous)
//            .withSubstances("UO2+2 UO2OH+ UO2(OH)2@ UO2CO3@ (UO2)3CO3(OH)3+ (UO2)2(OH)+3 UO2(CO3)3-4 UO2(CO3)2-2 "
//                            "(UO2)2(OH)2+2 (UO2)3(OH)5+ (UO2)4(OH)7+");
//    for(auto species : species_with_u)
//        std::cout << species.name() << std::endl;
//    getchar();

    // (UO2)2(OH)+3
    //(UO2)2(OH)2+2
    //(UO2)2CO3(OH)3-
    //(UO2)3(CO3)6-6
    //(UO2)3(OH)4+2
    //(UO2)3(OH)5+
    //(UO2)3(OH)7-
    //(UO2)3CO3(OH)3+
    //(UO2)4(OH)7+
    //Ca2UO2(CO3)3@
    //CaUO2(CO3)3-2
    //MgUO2(CO3)3-2
    //U(CO3)4-4
    //U(CO3)5-6
    //U(NO3)+3
    //U(NO3)2+2
    //U(OH)2+2
    //U(OH)3+
    //U(OH)4@
    //U(SCN)2+2
    //U(SO4)+2
    //U(SO4)2@
    //U+4
    //UCO3(OH)3-
    //UCl+3
    //UO2(CO3)2-2
    //UO2(CO3)3-4
    //UO2(CO3)3-5
    //UO2(H2PO4)+
    //UO2(H2PO4)2@
    //UO2(H3PO4)+2
    //UO2(HPO4)@
    //UO2(NO3)+
    //UO2(OH)2@
    //UO2(OH)3-
    //UO2(OH)4-2
    //UO2(PO4)-
    //UO2(SCN)2@
    //UO2(SCN)3-
    //UO2(SO4)2-2
    //UO2(SO4)3-4
    //UO2(SO4)@
    //UO2+
    //UO2+2
    //UO2CO3@
    //UO2Cl+
    //UO2Cl2@
    //UO2H5(PO4)2+
    //UO2HSiO3+
    //UO2OH+
    //UO2SCN+
    //UOH+3
    //USCN+3
    //"UO2+2 UO2OH+ UO2(OH)2@ UO2CO3@ (UO2)3CO3(OH)3+ (UO2)2(OH)+3 UO2(CO3)3-4"

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

    std::vector<double> pHs = {5, 6, 7, 8, 9, 10};
    std::vector<double> n(pHs);
    std::vector<double> nUO2(pHs);
    std::vector<double> nUO2OH(pHs);
    std::vector<double> nUO2OH2(pHs);
    std::vector<double> UO2CO3(pHs); // UO2CO3@
    std::vector<double> UO23CO3OH3(pHs); // (UO2)3CO3(OH)3+

    std::vector<double> UO22OH(pHs); // (UO2)2(OH)+3
    std::vector<double> UO2CO33(pHs); // UO2(CO3)3-4
    std::vector<double> UO2CO32(pHs); // UO2(CO3)2-2
    std::vector<double> UO22OH2(pHs); // (UO2)2(OH)2+2
    std::vector<double> UO23OH5(pHs); // (UO2)3(OH)5+
    std::vector<double> UO24OH7(pHs); // (UO2)4(OH)7+

    auto species_list = SpeciesList("UO2+2 UO2OH+ UO2(OH)2@ UO2CO3@ (UO2)3CO3(OH)3+ (UO2)2(OH)+3 "
                                    "UO2(CO3)3-4 UO2(CO3)2-2 (UO2)2(OH)2+2 (UO2)3(OH)5+ (UO2)4(OH)7+");
    auto bU = 0.0;

    // ---------------------------------------------------------------------------------------------------------------//
    // OrdinaryPortlandCement-PoreWater
    // ---------------------------------------------------------------------------------------------------------------//

    // Result file for the portland cement problem
    std::ofstream result_file;
    result_file.precision(6);
    result_file.open ("results-portlandcement.txt");

    // Define initial equilibrium state
    ChemicalState solutionstate_pc(system);
    solutionstate_pc.setSpeciesMass("H2O@",   1e3,  "g");
    solutionstate_pc.setSpeciesMass("CO2@",   1e-1, "g");
    solutionstate_pc.setSpeciesMass("SO2@",   1e-1, "g");
    solutionstate_pc.setSpeciesMass("H3PO4@", 1e-5, "g");

    // Aqueous properties of the chemical state
    AqueousProps aprops_pc(system);
    ChemicalProps props_pc(system);

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

    // Output the header of the table
    std::cout << "      pH       success       %n(UO2+2)   %n(UO2OH+) %n(UO2(OH)2@)  %n(UO2CO3@) "
                 "%n((UO2)3CO3(OH)3+) %n((UO2)2(OH)+3) %n(UO2(CO3)3-4) %n(UO2(CO3)2-2) %n((UO2)2(OH)2+2) %n((UO2)4(OH)7+)" << std::endl;

    for(Index i = 0; i < pHs.size(); i++)
    {
        conditions.pH(pHs[i]);

        // Equilibrate the initial state with given conditions and component amounts
        res = solver.solve(solutionstate_pc, conditions, portlandcementb);
        //std::cout << "res (OrdinaryPortlandCement-PoreWater) = " << res.optima.succeeded << std::endl;

        // Update aqueous properties
        aprops_pc.update(solutionstate_pc);
        props_pc.update(solutionstate_pc);
        //
        bU = props_pc.elementAmount("U");
        for(auto species : species_list)
        {
            result_file << solutionstate_pc.speciesAmount(species.name()) / bU * 100 << " ";
        }
        nUO2[i] = solutionstate_pc.speciesAmount("UO2+2");
        nUO2OH[i] = solutionstate_pc.speciesAmount("UO2OH+");
        nUO2OH2[i] = solutionstate_pc.speciesAmount("UO2(OH)2@");
        UO2CO3[i] = solutionstate_pc.speciesAmount("UO2CO3@");
        UO23CO3OH3[i] = solutionstate_pc.speciesAmount("(UO2)3CO3(OH)3+");
        UO22OH[i] = solutionstate_pc.speciesAmount("(UO2)2(OH)+3");
        UO2CO33[i] = solutionstate_pc.speciesAmount("UO2(CO3)3-4");
        UO2CO32[i] = solutionstate_pc.speciesAmount("UO2(CO3)2-2");
        UO22OH2[i] = solutionstate_pc.speciesAmount("(UO2)2(OH)2+2");
        UO24OH7[i] = solutionstate_pc.speciesAmount("(UO2)4(OH)7+");

        // Output results for the current pH
        std::cout.precision(6);
        std::cout << std::scientific << " "
                  << pHs[i] << "    "
                  << res.optima.succeeded << "        "
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
        solutionstate_pc.output("state-psinagra1207-portlandcement-pH-" + std::to_string(pHs[i]) + ".txt");
        result_file << "\n";
    }

    result_file.close();

    // ---------------------------------------------------------------------------------------------------------------//
    // Boda-Caly-pore-water
    // ---------------------------------------------------------------------------------------------------------------//

    // Result file for the portland cement problem
    std::ofstream result_file_bc;
    result_file_bc.precision(6);
    result_file_bc.open ("results-bodacaly.txt");

    // Define initial equilibrium state
    ChemicalState solutionstate_bc(system);
    solutionstate_bc.setSpeciesMass("H2O@",   1e3,  "g");
    solutionstate_bc.setSpeciesMass("CO2@",   1e-1, "g");
    solutionstate_bc.setSpeciesMass("SO2@",   1e-1, "g");
    solutionstate_bc.setSpeciesMass("H3PO4@", 1e-5, "g");

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
// Output the header of the table
    std::cout << "      pH       success       %n(UO2+2)   %n(UO2OH+) %n(UO2(OH)2@)  %n(UO2CO3@) "
                 "%n((UO2)3CO3(OH)3+) %n((UO2)2(OH)+3) %n(UO2(CO3)3-4) %n((UO2)2(OH)2+2) %n((UO2)4(OH)7+)" << std::endl;

    for(Index i = 0; i < pHs.size(); i++)
    {
        conditions.pH(pHs[i]);

        // Equilibrate the initial state with given conditions and component amounts
        res = solver.solve(solutionstate_bc, conditions, bodycalyb);
        //std::cout << "res (OrdinaryPortlandCement-PoreWater) = " << res.optima.succeeded << std::endl;

        // Update aqueous properties
        aprops_pc.update(solutionstate_bc);
        props_pc.update(solutionstate_bc);
        //
        bU = props_pc.elementAmount("U");
        for(auto species : species_list)
        {
            result_file << solutionstate_bc.speciesAmount(species.name()) / bU * 100 << " ";
        }
        nUO2[i] = solutionstate_bc.speciesAmount("UO2+2");
        nUO2OH[i] = solutionstate_bc.speciesAmount("UO2OH+");
        nUO2OH2[i] = solutionstate_bc.speciesAmount("UO2(OH)2@");
        UO2CO3[i] = solutionstate_bc.speciesAmount("UO2CO3@");
        UO23CO3OH3[i] = solutionstate_bc.speciesAmount("(UO2)3CO3(OH)3+");
        UO22OH[i] = solutionstate_bc.speciesAmount("(UO2)2(OH)+3");
        UO2CO33[i] = solutionstate_bc.speciesAmount("UO2(CO3)3-4");
        UO2CO32[i] = solutionstate_bc.speciesAmount("UO2(CO3)2-2");
        UO22OH2[i] = solutionstate_bc.speciesAmount("(UO2)2(OH)2+2");
        UO24OH7[i] = solutionstate_bc.speciesAmount("(UO2)4(OH)7+");

        // Output results for the current pH
        std::cout.precision(6);
        std::cout << std::scientific << " "
                  << pHs[i] << "    "
                  << res.optima.succeeded << "        "
                  << nUO2[i] / bU * 100 << " "
                  << nUO2OH[i] / bU * 100 << "  "
                  << nUO2OH2[i] / bU * 100 << " "
                  << UO2CO3[i] / bU * 100 << " "
                  << UO23CO3OH3[i] / bU * 100 << "       "
                  << UO22OH[i] / bU * 100 << " "
                  << UO2CO33[i] / bU * 100 << " "
                  << UO2CO32[i] / bU * 100 << " "
                  << UO22OH2[i] / bU * 100 << " "
                  << UO24OH7[i] / bU * 100 << std::endl;

        // Output the chemical state to a console
        solutionstate_pc.output("state-psinagra1207-bodycaly-pH-" + std::to_string(pHs[i]) + ".txt");
        result_file << "\n";
    }

    return 0;
}
