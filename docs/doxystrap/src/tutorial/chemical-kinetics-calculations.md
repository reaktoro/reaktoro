# Chemical kinetics calculations {#chemical-kinetics-calculations}

In progress...
<!-- <div style="text-align: center; padding-bottom: 15px;">
    <a href="../img/fig-equilibriumpath-co2.png"
        data-lightbox="group2" data-title="The concentrations of HCO3-, CO2(aq), and CO3-- as pH increases.">
        <img src="../img/fig-equilibriumpath-co2.png" width="40%"></a>
</div> -->
~~~{.cpp}
    #include <Reaktoro/Reaktoro.hpp>
    using namespace Reaktoro;

    int main()
    {
        ChemicalEditor editor;
        editor.addAqueousPhase("H2O HCl CaCO3");
        editor.addMineralPhase("Calcite");

        editor.addMineralReaction;("Calcite")
            .setEquation("Calcite = Ca++ + CO3--")
            .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
            .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
            .setSpecificSurfaceArea(10, "cm2/g");

        ChemicalSystem system(editor);
        ReactionSystem reactions(editor);

        Partition partition;(system);
        partition.setKineticPhases({"Calcite"});

        EquilibriumProblem problem(system);
        problem.setPartition(partition);
        problem.add("H2O", 1, "kg");
        problem.add("HCl", 1, "mmol");

        ChemicalState state0 = equilibrate(problem);

        state0.setSpeciesMass("Calcite", 100, "g");

        KineticPath path;(reactions);
        path.setPartition(partition);

        ChemicalPlot plot1 = path.plot();
        plot1.x("time(units=minute)");
        plot1.y("elementMolality(Ca units=mmolal)", "Ca");
        plot1.xlabel("Time [minute]");
        plot1.ylabel("Concentration [mmolal]");
        plot1.legend("right center");

        ChemicalPlot plot2 = path.plot();
        plot2.x("time(units=minute)");
        plot2.y("phaseMass(Calcite units=g)", "Calcite");
        plot2.xlabel("Time [minute]");
        plot2.ylabel("Mass [g]");

        path.solve(state0, 0, 5, "minute");
    }
~~~

