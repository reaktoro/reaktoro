// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

    ChemicalEditor editor(database);

    editor.addAqueousPhase({ "H2O(l)", "OH-", "H+", "HCO3-", "CO3--", "Fe++",
        "FeOH+", "FeOH++", "Fe+++", "H2(aq)", "O2(aq)", "HS-", "S5--", "S4--",
        "H2S(aq)", "S3--", "S2--", "SO4--", "HSO4-" }).setChemicalModelPitzerHMW();
    editor.addMineralPhase("Pyrrhotite");
    editor.addMineralPhase("Siderite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.setTemperature(25, "celsius");
    problem.setPressure(1, "atm");
    problem.add("H2O", 1.00, "mol");
    problem.add("H2S", 0.05, "mol");
    problem.add("Pyrrhotite", 0.0, "mol");
    problem.add("Siderite", 0.5, "mol");

    ChemicalState state = equilibrate(problem);

    std::cout << state << std::endl;
}
