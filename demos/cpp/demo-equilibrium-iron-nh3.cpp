// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
    editor.addAqueousPhaseWithElementsOf("H2O Fe(OH)2 Fe(OH)3 NH3");
    editor.addGaseousPhase("NH3(g)");
    editor.addMineralPhase("Magnetite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.add("H2O", 1, "kg");
    problem.add("Fe(OH)2", 1, "mol");
    problem.add("Fe(OH)3", 2, "mol");
    problem.add("NH3", 1, "mmol");

    ChemicalState state = equilibrate(problem);

    std::cout << state << std::endl;
}
