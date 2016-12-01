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
    editor.addAqueousPhase("H O C Ca Mg Na Cl");
    editor.addGaseousPhase("H2O(g) CO2(g)")
        .setChemicalModelSpycherPruessEnnis();
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Magnesite");
    editor.addMineralPhase("Dolomite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.setTemperature(60.0, "celsius");
    problem.setPressure(200.0, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CaCO3", 200, "g");
    problem.add("MgCO3", 50, "g");
    problem.add("NaCl", 0.1, "mol");

    EquilibriumState state1 = equilibrate(problem);

    problem.add("CO2", 2, "mol");

    EquilibriumState state2 = equilibrate(problem);

    EquilibriumPathOptions options;
    options.equilibrium.hessian = GibbsHessian::Exact;

    EquilibriumPath path(system);
    path.setOptions(options);

    ChemicalPlot plot0 = path.plot();
    plot0.x("t");
    plot0.y("Ca", "elementMolality(Ca)");
    plot0.y("Mg", "elementMolality(Mg)");
    plot0.title("Ca and Mg Concentration");
    plot0.xlabel("t");
    plot0.ylabel("Concentration [molal]");
    plot0.yformat("%g");
    plot0.legend("right center");

    ChemicalPlot plot1 = path.plot();
    plot1.x("t");
    plot1.y("pH", "pH");
    plot1.xlabel("t");
    plot1.ylabel("pH");

    ChemicalOutput output = path.output();
    output.data("t elementMolality(Ca) elementMolality(Mg) pH");
    output.headings("t Ca[molal] Mg[molal] ph");
    output.file("result.txt");

    path.solve(state1, state2);
}
