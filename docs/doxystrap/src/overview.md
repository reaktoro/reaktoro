# Overview {#mainpage}

Reaktoro is a framework developed in C++ and Python that implements numerical methods for modeling chemically reactive processes governed by either chemical equilibrium, chemical kinetics, or a combination of both.

## Chemical Equilibrium

Here is a simple C++ code using Reaktoro to perform a multiphase chemical equilibrium calculation:

~~~{.cpp}
    #include <Reaktoro/Reaktoro.hpp>
    using namespace Reaktoro;

    int main()
    {
        ChemicalEditor editor;
        editor.addAqueousPhase("H2O NaCl CaCO3 CO2");
        editor.addGaseousPhase("CO2(g)");
        editor.addMineralPhase("Calcite");

        ChemicalSystem system(editor);

        EquilibriumProblem problem(system);
        problem.add("H2O", 1, "kg");
        problem.add("CO2", 1, "mol");
        problem.add("NaCl", 0.7, "mol");
        problem.add("CaCO3", 1, "g");

        ChemicalState state = equilibrate(problem);

        state.output("result.txt");
    }
~~~

This calculation could also be performed using Reaktoro's Python interface:

~~~{.cpp}
using namespace Reaktoro; {delete}
from reaktoro import *

editor = ChemicalEditor()
; {delete}
ChemicalEditor editor; {delete}
editor.addAqueousPhase;("H2O NaCl CaCO3 CO2")
; {delete}
editor.addGaseousPhase;("CO2(g)")
; {delete}
editor.addMineralPhase;("Calcite")
; {delete}

ChemicalSystem system; {delete}
system = ChemicalSystem(editor)
; {delete}

EquilibriumProblem problem; {delete}
problem = EquilibriumProblem(system)
; {delete}
problem.add("H2O", 1, "kg")
; {delete}
problem.add("CO2", 1, "mol")
; {delete}
problem.add("NaCl", 0.7, "mol")
; {delete}
problem.add("CaCO3", 1, "g")
; {delete}

state = equilibrate(problem)

ChemicalState state; {delete}
state.output("result.txt")
~~~

@htmlonly
<a href="img/demo-equilibrium-mainpage.svg"
    data-fancybox="gallery-mainpage"
    data-caption="The output of the chemical equilibrium calculation.">
    <button type="button" class="btn btn-primary">Check the Result »</button>
</a>
@endhtmlonly

---

## Chemical Kinetics

%Reaktoro can also perform chemical kinetic calculations with both equilibrium-controlled and kinetically-controlled reactions. The C++ example below demonstrate this for a simple mineral dissolution modeling, in which CaCO3(s, calcite) reacts with a carbonated aqueous solution:

~~~{.cpp}
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
        plot1.y("speciesMass(Calcite units=g)", "Calcite");
        plot1.xlabel("Time [minute]");
        plot1.ylabel("Mass [g]");

        ChemicalPlot plot2 = path.plot();
        plot2.x("time(units=minute)");
        plot2.y("pH");
        plot2.xlabel("Time [minute]");
        plot2.ylabel("pH");

        ChemicalPlot plot3 = path.plot();
        plot3.x("time(units=minute)");
        plot3.y("speciesMolality(Ca++ units=mmolal)", "Ca++");
        plot3.y("speciesMolality(HCO3- units=mmolal)", "HCO3-");
        plot3.xlabel("Time [minute]");
        plot3.ylabel("Concentration [mmolal]");

        path.solve(initialstate, 0, 5, "minute");
    }
~~~

When the application is executed, the following figures are produced:

@htmlonly


<div class="container-fluid">
<div id="myCarousel" class="carousel slide" data-ride="carousel">
  <!-- Indicators -->
  <ol class="carousel-indicators" hidden>
    <li data-target="#myCarousel" data-slide-to="0" class="active"></li>
    <li data-target="#myCarousel" data-slide-to="1"></li>
    <li data-target="#myCarousel" data-slide-to="2"></li>
  </ol>

  <!-- Wrapper for slides -->
  <div class="carousel-inner">
    <div class="item active">
        <a href="img/demo-kineticpath-mainpage-1.svg"
            data-fancybox="gallery-kineticpath"
            data-caption="The mass of calcite mineral as it dissolves with time.">
            <img src="img/demo-kineticpath-mainpage-1.svg" style="height: 480px"></a>
    </div>

    <div class="item">
        <a href="img/demo-kineticpath-mainpage-2.svg"
            data-fancybox="gallery-kineticpath"
            data-caption="The pH of the aqueous phase as calcite dissolves.">
            <img src="img/demo-kineticpath-mainpage-2.svg" style="height: 480px"></a>
    </div>

    <div class="item">
        <a href="img/demo-kineticpath-mainpage-3.svg"
            data-fancybox="gallery-kineticpath"
            data-caption="The molalities of aqueous species Ca++ and HCO3- as calcite dissolves.">
            <img src="img/demo-kineticpath-mainpage-3.svg" style="height: 480px"></a>
    </div>
  </div>

  <!-- Left and right controls -->
  <a class="left carousel-control" href="#myCarousel" data-slide="prev">
    <span class="glyphicon glyphicon-chevron-left"></span>
    <span class="sr-only">Previous</span>
  </a>
  <a class="right carousel-control" href="#myCarousel" data-slide="next">
    <span class="glyphicon glyphicon-chevron-right"></span>
    <span class="sr-only">Next</span>
  </a>
</div>
</div>

@endhtmlonly

In the example above, the mineral reaction is specified to be under kinetic control and the aqueous species in chemical equilibrium at all times. As the mineral dissolves, it perturbs the chemical equilibrium state of the aqueous species. By assuming the aqueous species to be always in equilibrium, it is like if they were capable of reacting instantaneously to a new state of equilibrium. In general, the aqueous species react among themselves at much faster rates than mineral dissolution reactions, and thus this *partial equilibrium assumption* is plausible, and fairly accurate in most cases.

---

## Using %Reaktoro with GEMS and PHREEQC

There are excellent open-source software for modeling geochemical systems. Among all these software, two widely used are [GEMS] and [PHREEQC].

@htmlonly
<div class="container-fluid">
<div class="row">
<div class="col-sm-1 col-md-1"></div>
<div class="col-sm-4 col-md-4">
<a href="http://gems.web.psi.ch/" target="_blank">
<img src="img/gems-logo.png" width="100%"/></a>
</div>
<div class="col-sm-1 col-md-1"></div>
<div class="col-sm-1 col-md-1"></div>
<div class="col-sm-4 col-md-4">
<a href="https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/" target="_blank"> <img src="img/phreeqc-logo.svg" width="100%"/></a>
</div>
<div class="col-md-1"></div>
</div>
</div>
@endhtmlonly

Because many users rely on specific thermodynamic models (e.g., models for activity coefficients, fugacity coefficients) and databases (e.g., PSI/NAGRA, PHREEQC databases) supported in these two software, two classes have been implemented in Reaktoro to permit the combined use of those models and databases with Reaktoro's numerical methods. These classes are [Gems](@ref Reaktoro::Gems) and [Phreeqc](@ref Reaktoro::Phreeqc).

<!-- In what follows, a [ChemicalSystem](@ref Reaktoro::ChemicalSystem) object is constructed by converting
[Gems](@ref Reaktoro::Gems) and [Phreeqc](@ref Reaktoro::Phreeqc)
a CBy using These will be solved using Reaktoro's numerical methods.
GEMS will be used to calculate thermodynamic properties
of species (e.g., activities, standard chemical potentials) and
phases (e.g., molar volume, heat capacity, etc.). -->

The code below briefly demonstrates how GEMS thermodynamic models and databases can be used with Reaktoro's numerical methods.

~~~{.cpp}
    #include <Reaktoro/Reaktoro.hpp>
    using namespace Reaktoro;

    int main()
    {
        // Use an exported project file from GEMS to initialize a Gems object,
        Gems gems;("project.lst");

        // and then use it to construct the ChemicalSystem object.
        ChemicalSystem system = gems;

        // Create a ChemicalState object that contains the temperature, pressure,
        // and amounts of species stored in the exported GEMS file.
        ChemicalState state = gems.state(system);

        // Change the temperature of the chemical state,
        state.setTemperature(50, "celsius");

        // and then equilibrate the modified chemical state using Reaktoro's methods.
        equilibrate(state);

        // Output the updated equilibrium state to a file.
        state.output("state.txt");
    }
~~~

Similarly, the code below demonstrate the combined use of Reaktoro and [PHREEQC].

~~~{.cpp}
    #include <Reaktoro/Reaktoro.hpp>
    using namespace Reaktoro;

    int main()
    {
        // Initialize a Phreeqc object with a PHREEQC database file such as
        // phreeqc.dat, pitzer.dat, llnl.dat, sit.dat, minteq.dat, etc.
        Phreeqc phreeqc;("some-phreeqc-database-file.dat");

        // Use the just created Phreeqc object to execute a PHREEQC script,
        phreeqc.execute("some-phreeqc-input-script.dat");

        // and then use it to construct a ChemicalSystem object containing all
        // phases and species selected by PHREEQC.
        ChemicalSystem system = phreeqc;

        // Create a ChemicalState object for the PHREEQC calculated chemical state,
        // which was produced when the Phreeqc object executed the script above.
        ChemicalState state = phreeqc.state(system);

        // Define an equilibrium problem that uses the calculated PHREEQC state as
        // a basis for the calculation of a new equilibrium state. The problem
        // below adds 0.1 moles of CO2 to the calculated PHREEQC state while
        // maintaining the same temperature and pressure specified in the script.
        EquilibriumProblem problem(system);
        problem.setTemperature(state.temperature());
        problem.setPressure(state.pressure());
        problem.add(state);
        problem.add("CO2", 0.1, "moles");

        // Calculate the new equilibrium state using the above equilibrium problem.
        equilibrate(state, problem);

        // Output the updated equilibrium state to a file.
        state.output("state.txt");
    }
~~~

--------------------------------
<!-- 
Thermodynamic Property Calculations
-----------------------------------

~~~{.cpp}

~~~ -->

[PHREEQC]: http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
[GEMS]: http://gems.web.psi.ch/