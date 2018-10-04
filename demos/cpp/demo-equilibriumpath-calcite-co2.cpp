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
    ChemicalEditor editor;
    editor.addAqueousPhase("H O C Ca Mg Na Cl");
    editor.addGaseousPhase({"H2O(g)", "CO2(g)"})
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

    ChemicalState state1 = equilibrate(problem);

    problem.add("CO2", 2, "mol");

    ChemicalState state2 = equilibrate(problem);

    EquilibriumPathOptions options;
    options.equilibrium.hessian = GibbsHessian::Exact;

    EquilibriumPath path(system);
    path.setOptions(options);

    ChemicalPlot plot0 = path.plot();
    plot0.x("t");
    plot0.y("elementMolality(Ca)", "Ca");
    plot0.y("elementMolality(Mg)", "Mg");
    plot0.title("Ca and Mg Concentration");
    plot0.xlabel("t");
    plot0.ylabel("Concentration [molal]");
    plot0.yformat("%g");
    plot0.legend("right center");

    ChemicalPlot plot1 = path.plot();
    plot1.x("t");
    plot1.y("pH");
    plot1.xlabel("t");
    plot1.ylabel("pH");

    ChemicalOutput output = path.output();
    output.filename("result.txt");
    output.add("t");
    output.add("elementMolality(Ca)", "Ca [molal]");
    output.add("elementMolality(Mg)", "Mg [molal]");
    output.add("pH");

    path.solve(state1, state2);
}
