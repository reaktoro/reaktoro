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
    ChemicalEditor editor;
    editor.addAqueousPhase("H O C Ca Cl Mg");
    editor.addGaseousPhase("H2O(g) CO2(g) H2(g) O2(g) CH4(g)");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Dolomite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.add("H2O", 1, "kg");
    problem.add("MgCl2", 0.002, "mol");
    problem.add("CaCO3", 1, "mol");

    EquilibriumState state = equilibrate(problem);

    std::cout << state << std::endl;
}
