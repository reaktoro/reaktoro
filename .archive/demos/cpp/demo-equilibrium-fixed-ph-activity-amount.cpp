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
    ChemicalEditor editor;
    editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C");
    editor.addGaseousPhaseWithElements("H O C");

    ChemicalSystem system(editor);

    EquilibriumInverseProblem problem(system);
    problem.add("H2O", 1, "kg");
    problem.add("NaCl", 0.1, "mol");
    problem.add("CaCl2", 2, "mmol");
    problem.add("MgCl2", 4, "mmol");
    problem.pH(3.0, "HCl");
    problem.fixSpeciesAmount("CO2(g)", 1.0, "mol");
    problem.fixSpeciesActivity("O2(g)", 0.20);

    ChemicalState state = equilibrate(problem);

    std::cout << state << std::endl;
}
