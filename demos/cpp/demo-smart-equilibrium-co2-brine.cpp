// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

int main()
{
    Database database("supcrt98.xml");

    ChemicalEditor editor(database);
    editor.addAqueousPhaseWithElements("H O Na Cl C Ca Mg Si")
        // .setChemicalModelPitzerHMW()
        // .setChemicalModelDebyeHuckel()
        // .setActivityModelDrummondCO2()
        ;
    editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
    editor.addMineralPhase("Halite");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Magnesite");
    editor.addMineralPhase("Aragonite");
    editor.addMineralPhase("Dolomite");
    editor.addMineralPhase("Quartz");

    ChemicalSystem system(editor);

    Partition partition(system);
    // partition.setInertSpecies({{"Halite"}, {"Quartz"}});

    EquilibriumProblem problem(partition);
    problem.setTemperature(60, "celsius");
    problem.setPressure(300, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CO2", 2, "mol");
    problem.add("NaCl", 0.1, "mol");
    // problem.add("MgCl2", 0.1, "mol");
    problem.add("Calcite", 1, "mol");
    // problem.add("Dolomite", 1, "mol");
    problem.add("Quartz", 1, "mol");

    EquilibriumResult result;
    SmartEquilibriumResult smartresult;

    EquilibriumSolver solver(partition);
    SmartEquilibriumSolver smartsolver(partition);

    std::cout << "Solving Problem 1...";
    ChemicalState state1(system);
    ChemicalState smartstate1(system);
    result = solver.solve(state1, problem);
    smartresult = smartsolver.solve(smartstate1, problem);
    std::cout << "Smart Prediction Succeess? " << (smartresult.estimate.accepted ? "yes" : "no") << std::endl;

    problem.add("MgCl2", 0.01, "mol");
    std::cout << "Solving Problem 2...";
    ChemicalState state2(system);
    ChemicalState smartstate2(system);
    result = solver.solve(state2, problem);
    smartresult = smartsolver.solve(smartstate2, problem);
    std::cout << "Smart Prediction Succeess? " << (smartresult.estimate.accepted ? "yes" : "no") << std::endl;

    problem.add("MgCl2", 0.01, "mol");
    std::cout << "Solving Problem 3...";
    ChemicalState state3(system);
    ChemicalState smartstate3(system);
    result = solver.solve(state3, problem);
    smartresult = smartsolver.solve(smartstate3, problem);
    std::cout << "Smart Prediction Succeess? " << (smartresult.estimate.accepted ? "yes" : "no") << std::endl;

    problem.add("MgCl2", 0.01, "mol");
    std::cout << "Solving Problem 4...";
    ChemicalState state4(system);
    ChemicalState smartstate4(system);
    result = solver.solve(state4, problem);
    smartresult = smartsolver.solve(smartstate4, problem);
    std::cout << "Smart Prediction Succeess? " << (smartresult.estimate.accepted ? "yes" : "no") << std::endl;

    state1.output("state1.txt");
    state2.output("state2.txt");
    state3.output("state3.txt");
    state4.output("state4.txt");

    smartstate1.output("smartstate1.txt");
    smartstate2.output("smartstate2.txt");
    smartstate3.output("smartstate3.txt");
    smartstate4.output("smartstate4.txt");

    const Vector n1 = state1.speciesAmounts();
    const Vector n2 = state2.speciesAmounts();
    const Vector n3 = state3.speciesAmounts();
    const Vector n4 = state4.speciesAmounts();

    const Vector y1 = state1.equilibrium().elementChemicalPotentials();
    const Vector y2 = state2.equilibrium().elementChemicalPotentials();
    const Vector y3 = state3.equilibrium().elementChemicalPotentials();
    const Vector y4 = state4.equilibrium().elementChemicalPotentials();

    const Vector z1 = state1.equilibrium().speciesStabilities();
    const Vector z2 = state2.equilibrium().speciesStabilities();
    const Vector z3 = state3.equilibrium().speciesStabilities();
    const Vector z4 = state4.equilibrium().speciesStabilities();

    const Vector smart_n1 = smartstate1.speciesAmounts();
    const Vector smart_n2 = smartstate2.speciesAmounts();
    const Vector smart_n3 = smartstate3.speciesAmounts();
    const Vector smart_n4 = smartstate4.speciesAmounts();

    const Vector smart_y1 = smartstate1.equilibrium().elementChemicalPotentials();
    const Vector smart_y2 = smartstate2.equilibrium().elementChemicalPotentials();
    const Vector smart_y3 = smartstate3.equilibrium().elementChemicalPotentials();
    const Vector smart_y4 = smartstate4.equilibrium().elementChemicalPotentials();

    const Vector smart_z1 = smartstate1.equilibrium().speciesStabilities();
    const Vector smart_z2 = smartstate2.equilibrium().speciesStabilities();
    const Vector smart_z3 = smartstate3.equilibrium().speciesStabilities();
    const Vector smart_z4 = smartstate4.equilibrium().speciesStabilities();

    const Vector ndiff1 = abs(n1 - smart_n1);
    const Vector ndiff2 = abs(n2 - smart_n2);
    const Vector ndiff3 = abs(n3 - smart_n3);
    const Vector ndiff4 = abs(n4 - smart_n4);

    const Vector ydiff1 = abs(y1 - smart_y1);
    const Vector ydiff2 = abs(y2 - smart_y2);
    const Vector ydiff3 = abs(y3 - smart_y3);
    const Vector ydiff4 = abs(y4 - smart_y4);

    const Vector zdiff1 = abs(z1 - smart_z1);
    const Vector zdiff2 = abs(z2 - smart_z2);
    const Vector zdiff3 = abs(z3 - smart_z3);
    const Vector zdiff4 = abs(z4 - smart_z4);

    Index ispecies1;
    Index ispecies2;
    Index ispecies3;
    Index ispecies4;

    const double max_ndiff1 = ndiff1.maxCoeff(&ispecies1);
    const double max_ndiff2 = ndiff2.maxCoeff(&ispecies2);
    const double max_ndiff3 = ndiff3.maxCoeff(&ispecies3);
    const double max_ndiff4 = ndiff4.maxCoeff(&ispecies4);

    const double max_ydiff1 = ydiff1.maxCoeff();
    const double max_ydiff2 = ydiff2.maxCoeff();
    const double max_ydiff3 = ydiff3.maxCoeff();
    const double max_ydiff4 = ydiff4.maxCoeff();

    const double max_zdiff1 = zdiff1.maxCoeff();
    const double max_zdiff2 = zdiff2.maxCoeff();
    const double max_zdiff3 = zdiff3.maxCoeff();
    const double max_zdiff4 = zdiff4.maxCoeff();

    std::cout << std::scientific;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(25) << "Species";
    std::cout << std::left << std::setw(25) << "|smart(n1) - conv(n1)|";
    std::cout << std::left << std::setw(25) << "|smart(n2) - conv(n2)|";
    std::cout << std::left << std::setw(25) << "|smart(n3) - conv(n3)|";
    std::cout << std::left << std::setw(25) << "|smart(n4) - conv(n4)|";
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    for(auto i = 0; i < system.numSpecies(); ++i)
    {
        std::cout << std::left << std::setw(25) << system.species(i).name();
        std::cout << std::left << std::setw(25) << ndiff1[i];
        std::cout << std::left << std::setw(25) << ndiff2[i];
        std::cout << std::left << std::setw(25) << ndiff3[i];
        std::cout << std::left << std::setw(25) << ndiff4[i];
        std::cout << std::endl;
    }
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(25) << "Element";
    std::cout << std::left << std::setw(25) << "|smart(y1) - conv(y1)|";
    std::cout << std::left << std::setw(25) << "|smart(y2) - conv(y2)|";
    std::cout << std::left << std::setw(25) << "|smart(y3) - conv(y3)|";
    std::cout << std::left << std::setw(25) << "|smart(y4) - conv(y4)|";
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    for(auto i = 0; i < system.numElements(); ++i)
    {
        std::cout << std::left << std::setw(25) << system.element(i).name();
        std::cout << std::left << std::setw(25) << ydiff1[i] / (universalGasConstant * state1.temperature());
        std::cout << std::left << std::setw(25) << ydiff2[i] / (universalGasConstant * state2.temperature());
        std::cout << std::left << std::setw(25) << ydiff3[i] / (universalGasConstant * state3.temperature());
        std::cout << std::left << std::setw(25) << ydiff4[i] / (universalGasConstant * state4.temperature());
        std::cout << std::endl;
    }
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(25) << "Species";
    std::cout << std::left << std::setw(25) << "|smart(z1) - conv(z1)|";
    std::cout << std::left << std::setw(25) << "|smart(z2) - conv(z2)|";
    std::cout << std::left << std::setw(25) << "|smart(z3) - conv(z3)|";
    std::cout << std::left << std::setw(25) << "|smart(z4) - conv(z4)|";
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    for(auto i = 0; i < system.numSpecies(); ++i)
    {
        std::cout << std::left << std::setw(25) << system.species(i).name();
        std::cout << std::left << std::setw(25) << zdiff1[i] / (universalGasConstant * state1.temperature());
        std::cout << std::left << std::setw(25) << zdiff2[i] / (universalGasConstant * state2.temperature());
        std::cout << std::left << std::setw(25) << zdiff3[i] / (universalGasConstant * state3.temperature());
        std::cout << std::left << std::setw(25) << zdiff4[i] / (universalGasConstant * state4.temperature());
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Max deviation of conventional and smart solutions of problem 1 = " << max_ndiff1/n1[ispecies1] * 100 << "% by species " << system.species(ispecies1).name() << std::endl;
    std::cout << "Max deviation of conventional and smart solutions of problem 2 = " << max_ndiff2/n2[ispecies2] * 100 << "% by species " << system.species(ispecies2).name() << std::endl;
    std::cout << "Max deviation of conventional and smart solutions of problem 3 = " << max_ndiff3/n3[ispecies3] * 100 << "% by species " << system.species(ispecies3).name() << std::endl;
    std::cout << "Max deviation of conventional and smart solutions of problem 4 = " << max_ndiff4/n4[ispecies4] * 100 << "% by species " << system.species(ispecies4).name() << std::endl;
}
