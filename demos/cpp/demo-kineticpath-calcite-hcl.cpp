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
    editor.addAqueousPhase("H2O HCl CaCO3");
    editor.addMineralPhase("Calcite");
//
    editor.addMineralReaction("Calcite")
        .setEquation("Calcite = Ca++ + CO3--")
        .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
        .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
        .setSpecificSurfaceArea(10, "cm2/g");

    ChemicalSystem system(editor);
    ReactionSystem reactions(editor);

    Partition partition(system);
//    partition.setKineticPhases({"Calcite"});

    EquilibriumProblem problem(system);
    problem.setPartition(partition);
    problem.add("H2O", 1, "kg");
    problem.add("HCl", 1, "mmol");
    problem.add("CaCO3", 1, "mmol");

    KineticState state0 = equilibrate(problem);

    state0.scalePhaseVolume("Aqueous", 1.0);

//    state0.setSpeciesMass("Calcite", 100, "g");
//
//    state0.scalePhaseVolume("Aqueous", 0.9);
//    state0.scalePhaseVolume("Calcite", 0.1);

    KineticOptions options;
    options.equilibrium.hessian = GibbsHessian::Exact;
    options.equilibrium.optimum.output = true;

    KineticPath path(reactions);
    path.setOptions(options);
    path.setPartition(partition);
    path.addFluidSink(0.1, "m3/s");

//    ChemicalPlot plot1 = path.plot();
//    plot1.x("time(units=minute)");
//    plot1.y("elementMolality(Ca)");
//    plot1.legend("Ca");
//    plot1.xlabel("Time [minute]");
//    plot1.ylabel("Concentration [molal]");
//    plot1.key("right center");
//
//    ChemicalPlot plot2 = path.plot();
//    plot2.x("time(units=minute)");
//    plot2.y("phaseMass(Calcite units=g)");
//    plot2.legend("Calcite");
//    plot2.xlabel("Time [minute]");
//    plot2.ylabel("Mass [g]");
//
//    ChemicalPlot plot3 = path.plot();
//    plot3.x("time(units=second)");
//    plot3.y("fluidVolume");
//    plot3.nolegend();
//    plot3.xlabel("Time [second]");
//    plot3.ylabel("Fluid Volume");
//
//    ChemicalPlot plot3 = path.plot();
//    plot3.x("time(units=minute)");
//    plot3.y("porosity");
//    plot3.nolegend();
//    plot3.xlabel("Time [minute]");
//    plot3.ylabel("Porosity");

    ChemicalOutput output = path.output();
    output.terminal(true);
    output.data("t elementMolality(Ca) phaseMass(Calcite)");

    path.solve(state0, 0, 5, "minute");
}
