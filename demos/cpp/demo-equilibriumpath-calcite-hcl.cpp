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
    editor.addAqueousPhase("H O Ca C Cl");
    editor.addMineralPhase("Calcite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem1(system);
    problem1.setTemperature(30.0, "celsius");
    problem1.setPressure(1.0, "bar");
    problem1.add("H2O", 1, "kg");
    problem1.add("CaCO3", 1, "g");

    EquilibriumProblem problem2(system);
    problem2.setTemperature(30.0, "celsius");
    problem2.setPressure(1.0, "bar");
    problem2.add("H2O", 1, "kg");
    problem2.add("CaCO3", 1, "g");
    problem2.add("HCl", 1, "mmol");

    ChemicalState state1 = equilibrate(problem1);
    ChemicalState state2 = equilibrate(problem2);

    EquilibriumPath path(system);

    ChemicalPlot plot1 = path.plot();
    plot1.x("elementAmount(Cl units=mmol)");
    plot1.y("pH");
    plot1.xlabel("HCl [mmol]");
    plot1.ylabel("pH");
    plot1.showlegend(false);

    ChemicalPlot plot2 = path.plot();
    plot2.x("elementAmount(Cl units=mmol)");
    plot2.y("Ca", "elementMolality(Ca units=mmolal)");
    plot2.xlabel("HCl [mmol]");
    plot2.ylabel("Concentration [mmolal]");
    plot2.legend("right center");

    ChemicalPlot plot3 = path.plot();
    plot3.x("elementAmount(Cl units=mmol)");
    plot3.y("CO2(aq)", "speciesMolality(CO2(aq) units=mmolal)");
    plot3.y("CO3--", "speciesMolality(CO3-- units=mmolal)");
    plot3.xlabel("HCl [mmol]");
    plot3.ylabel("Concentration [mmolal]");
    plot3.legend("right top");

    ChemicalPlot plot4 = path.plot();
    plot4.x("elementAmount(Cl units=mmol)");
    plot4.y("Calcite", "speciesMass(Calcite units=g)");
    plot4.xlabel("HCl [mmol]");
    plot4.ylabel("Mass [g]");

    ChemicalOutput output = path.output();
    output.file("result.txt");
    output.add("Cl [mmol]", "elementAmount(Cl units=mmol)");
    output.add("Ca [mmolal]", "elementMolality(Ca units=mmolal)");
    output.add("pH");
    output.add("speciesMass(Calcite units=g)");

    path.solve(state1, state2);
}
