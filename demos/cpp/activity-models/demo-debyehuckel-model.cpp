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
    DebyeHuckelParams db;
    db.setPHREEQC();

    ChemicalEditor editor;
    editor.addAqueousPhase("H2O NaCl CaCl2 CO2")
        .setChemicalModelDebyeHuckel(db);

    ChemicalSystem system(editor);

    EquilibriumProblem problem1(system);
    problem1.add("H2O", 1, "kg");

    EquilibriumProblem problem2(system);
    problem2.add("H2O", 1, "kg");
    problem2.add("NaCl", 0.1, "mol");
    problem2.add("CaCl2", 0.5, "mol");
    problem2.add("CO2", 0.2, "mol");

    EquilibriumState state1 = equilibrate(problem1);
    EquilibriumState state2 = equilibrate(problem2);

    EquilibriumPath path(system);

    ChemicalPlot plot = path.plot();
    plot.x("ionicStrength");
    plot.y("Na+",  "activityCoefficient(Na+)");
    plot.y("Cl-",  "activityCoefficient(Cl-)");
    plot.y("Ca++", "activityCoefficient(Ca++)");
    plot.ylabel("Activity Coefficient");
    plot.xlabel("I [molal]");
    path.solve(state1, state2);
}
