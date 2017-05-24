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
    editor.addAqueousPhase("H2O CO2 CaCO3");
    editor.addMineralPhase("Calcite");

    MineralReaction reaction = editor.addMineralReaction("Calcite");
    reaction.setEquation("Calcite = Ca++ + CO3--");
    reaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol");
    reaction.addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0");
    reaction.setSpecificSurfaceArea(10, "cm2/g");

    ChemicalSystem system(editor);
    ReactionSystem reactions(editor);

    Partition partition(system);
    partition.setKineticSpecies({"Calcite"});

    EquilibriumProblem problem(system);
    problem.setPartition(partition);
    problem.setTemperature(60, "celsius");
    problem.setPressure(100, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CO2", 0.1, "mol");

    ChemicalState initialstate = equilibrate(problem);

    initialstate.setSpeciesMass("Calcite", 100, "g");

    KineticPath path(reactions);
    path.setPartition(partition);

    ChemicalPlot plot1 = path.plot();
    plot1.x("time(units=minute)");
    plot1.y("Calcite", "speciesMass(Calcite units=g)");
    plot1.xlabel("Time [minute]");
    plot1.ylabel("Mass [g]");

    ChemicalPlot plot2 = path.plot();
    plot2.x("time(units=minute)");
    plot2.y("pH");
    plot2.xlabel("Time [minute]");
    plot2.ylabel("pH");

    ChemicalPlot plot3 = path.plot();
    plot3.x("time(units=minute)");
    plot3.y("Ca++", "speciesMolality(Ca++ units=mmolal)");
    plot3.y("HCO3-", "speciesMolality(HCO3- units=mmolal)");
    plot3.xlabel("Time [minute]");
    plot3.ylabel("Concentration [mmolal]");

    path.solve(initialstate, 0, 5, "minute");
}
