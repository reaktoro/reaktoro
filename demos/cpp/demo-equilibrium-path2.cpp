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

    EquilibriumOptions options;
    options.optimum.output.active = true;

    ChemicalState state1 = equilibrate(problem);

    problem.add("CO2", 2, "mol");

    ChemicalState state2 = equilibrate(problem);

    EquilibriumPathOptions pathoptions;
    pathoptions.equilibrium.hessian = GibbsHessian::Exact;

    EquilibriumPath path(system);
    path.setOptions(pathoptions);

    auto plots = path.plots(2);

    plots[0].xdata("t");
    plots[0].ydata("molality element=Ca units=molal;"
               "molality element=Mg units=molal;");
    plots[0] << "set xlabel 't'";
    plots[0] << "set ylabel 'Concentration [molal]'";
    plots[0] << "set format x '%g'";

    plots[1].xdata("t");
    plots[1].ydata("pH");
    plots[1] << "set xlabel 't'";
    plots[1] << "set ylabel 'pH'";

    ChemicalOutput output = path.output();
    output.header("t; Ca [molal]; Mg [molal]; pH");
    output.data("t; molality element=Ca; molality element=Mg; pH");
    output.file("result.txt");

    path.solve(state1, state2);
}




