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
    editor.addAqueousPhase("H O C Na Cl");

    ChemicalSystem system(editor);

    EquilibriumProblem problem1(system);
    problem1.add("H2O", 1, "kg");
    problem1.add("CO2", 0.5, "mol");
    problem1.add("HCl", 1, "mol");

    EquilibriumProblem problem2(system);
    problem2.add("H2O", 1, "kg");
    problem2.add("CO2", 0.5, "mol");
    problem2.add("NaOH", 2, "mol");

    ChemicalState state1 = equilibrate(problem1);
    ChemicalState state2 = equilibrate(problem2);

    EquilibriumPath path(system);

    ChemicalPlot plot = path.plot();
    plot.x("pH");
    plot.y("speciesMolality(HCO3-) speciesMolality(CO2(aq)) speciesMolality(CO3--)");
    plot.legend("HCO@_3^- CO_2(aq) CO@_3^{2-}");
    plot.xlabel("pH");
    plot.ylabel("Concentration [molal]");
    plot.yformat("%g");
    plot.key("left center Left reverse");

    ChemicalOutput output = path.output();
    output.data("t pH speciesMolality(HCO3-) speciesMolality(CO2(aq)) speciesMolality(CO3--)");
    output.headings("t pH HCO3-[molal] CO2(aq)[molal] CO3--[molal]");
    output.file("result.txt");

    path.solve(state1, state2);
}
