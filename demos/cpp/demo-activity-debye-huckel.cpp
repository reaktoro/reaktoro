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

//TODO - fix convergence error
int main()
{
    DebyeHuckelParams db;
    db.setPHREEQC();

    ChemicalEditor editor;
    editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCl2 CO2")
        .setChemicalModelDebyeHuckel(db);

    ChemicalSystem system(editor);

    EquilibriumProblem problem1(system);
    problem1.add("H2O", 1, "kg");

    EquilibriumProblem problem2(system);
    problem2.add("H2O", 1, "kg");
    problem2.add("NaCl", 0.1, "mol");
    problem2.add("CaCl2", 0.5, "mol");
    problem2.add("CO2", 0.2, "mol");

    ChemicalState state1 = equilibrate(problem1);
    ChemicalState state2 = equilibrate(problem2);

    EquilibriumPath path(system);

    ChemicalPlot plot = path.plot();
    plot.x("ionicStrength");
    plot.y("activityCoefficient(Na+)", "Na+");
    plot.y("activityCoefficient(Cl-)", "Cl-");
    plot.y("activityCoefficient(Ca++)", "Ca++");
    plot.ylabel("Activity Coefficient");
    plot.xlabel("I [molal]");
    path.solve(state1, state2);
}
