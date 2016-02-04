// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database database("databases/supcrt/supcrt98.xml");

    ChemicalEditor editor(database);
    editor.addAqueousPhase("H2O NaCl CaCO3");
    editor.addGaseousPhase("H2O(g) CO2(g)");
    editor.addMineralPhase("Calcite");

    ChemicalSystem system(editor);

    std::cout << system << std::endl;

    EquilibriumProblem problem(system);
    problem.add("H2O", 1, "kg");
    problem.add("NaCl", 0.1, "mol");
    problem.add("CaCO3", 10.0, "mol");
//    problem.add("HCl", 2.0, "mol");
//    problem.add("CO2", 10, "moles");
//    problem.setSpeciesAmount("Calcite", 10, "moles");
//    problem.setSpeciesAmount("CO2(g)", 10, "moles");
//    problem.setSpeciesAmount("H2O(l)", 55, "moles");
//    problem.setPhaseAmount("Aqueous", 60, "moles", "H2O");
//    problem.setPhaseAmount("Aqueous", 60, "moles");
//    problem.setPhaseAmount("Calcite", 10, "moles");
    problem.pH(3.0, "HCl");
//    problem.pH(11.0, "HCl", "NaOH");
//    problem.setPhaseVolume("Aqueous", 0.5, "m3");
//    problem.setPhaseVolume("Calcite", 0.5, "m3");

    EquilibriumOptions options;
//    options.optimum.output.active = true;
    options.hessian = GibbsHessian::Exact;

//    ChemicalState state = equilibrate(problem, "output=true");
    ChemicalState state = equilibrate(problem, options);

    std::cout << state << std::endl;
}
