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
    
    

    Database database("supcrt07.xml");
    
    DebyeHuckelParams dhModel{};
    //dhModel.setPHREEQC();
    dhModel.setKielland1937();

    ChemicalEditor editor(database);
    
    //editor.addAqueousPhase("H2O(l) FeS FeS2 H2S FeCO3").setChemicalModelDebyeHuckel(db);
    //editor.addAqueousPhase({ "H2O(l)","H+","OH-", "H2(aq)", "O2(aq)", "HS-", "H2S(aq)", "SO4--", "HSO4-" });
    //editor.addAqueousPhase("H2O H2S FeCO3");
    //editor.addAqueousPhase("H2O(l)");
    //editor.addMineralPhase("Pyrite");
    //editor.addMineralPhase("Pyrrhotite");
    //editor.addMineralPhase("Siderite");
    
    //editor.addAqueousPhase("H2O(l) FeS2 FeS H2S FeCO3 FeCO3(aq)").setChemicalModelDebyeHuckel(dhModel);
    //editor.addAqueousPhase("H2O(l) FeS2 FeS H2S FeCO3 FeCO3(aq)").setChemicalModelDebyeHuckel(dhModel);
    //editor.addAqueousPhase({ "H2O(l)", "OH-", "H+", "HCO3-", "CO3--", "Fe++", "FeOH+", "FeOH++", "Fe+++", "H2(aq)", "O2(aq)", "HS-", "S5--", "S4--", "H2S(aq)", "S3--", "S2--", "SO4--", "HSO4-" }).setChemicalModelDebyeHuckel(dhModel);
    editor.addAqueousPhase({ "H2O(l)", "OH-", "H+", "HCO3-", "CO3--", "Fe++", "FeOH+", "FeOH++", "Fe+++", "H2(aq)", "O2(aq)", "HS-", "S5--", "S4--", "H2S(aq)", "S3--", "S2--", "SO4--", "HSO4-" }).setChemicalModelPitzerHMW();
    //editor.addAqueousPhase("H2O(l) FeS2 FeS H2S FeCO3 FeCO3(aq)").setChemicalModelHKF();
    //editor.addAqueousPhase("H2O(l) FeS2 FeS H2S FeCO3").setChemicalModelPitzerHMW();
    //editor.addMineralPhase("Pyrite");
    editor.addMineralPhase("Pyrrhotite");
    editor.addMineralPhase("Siderite");


    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.setTemperature(25, "celsius");
    problem.setPressure(1, "atm");
    problem.add("H2O", 58.0, "kg");
    problem.add("HS-", 0.05, "mol");
    //problem.add("Pyrite", 0.0, "mol");
    problem.add("Pyrrhotite", 0.0, "mol");
    problem.add("Siderite", 0.5, "mol");
    //problem.pH(5.728);
    //problem.pE(-2.518);
    

    ChemicalState state = equilibrate(problem);

    /*
    problem.add("Siderite", 0.1, "mol");

    equilibrate(state, problem);
    */

    state.output("siderite_out.txt");

    std::cout << state << std::endl;
    std::cout << "RUN" << std::endl;
}
