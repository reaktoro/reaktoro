@mainpage notitle

Quick Start
===========

Reaktoro is a framework developed in C++ and Python that implements numerical methods for modeling chemically reactive processes governed by either chemical equilibrium, chemical kinetics, or a combination of both.

Example: Chemical Equilibrium
------------------------------
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

    EquilibriumState state = equilibrate(problem);

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

EquilibriumState state; {delete}
state.output("result.txt")
~~~

@htmlonly
<a href="img/demo-equilibrium-mainpage.svg"
    data-fancybox="gallery-mainpage"
    data-caption="The output of the chemical equilibrium calculation.">
    <button type="button" class="btn btn-primary">Check the Result »</button>
</a>
@endhtmlonly

--------------------------------

Example: Chemical Kinetics
--------------------------
Reaktoro can also perform chemical kinetic calculations with both equilibrium-controlled and kinetically-controlled reactions:

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

    KineticState initialstate = equilibrate(problem);

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

------------------------

Tutorials
---------

The two examples above show two major numerical capabilities of Reaktoro: *chemical equilibrium and kinetic calculations*. To learn more about how to use Reaktoro for these calculations, we have prepared [tutorials](@ref tutorials) to walk you through. These tutorials will explain step-by-step the use of each component in Reaktoro, which include class such as [ChemicalSystem](@ref Reaktoro::ChemicalSystem), [EquilibriumProblem](@ref Reaktoro::EquilibriumProblem), and [KineticPath](@ref Reaktoro::KineticPath), and functions, such as the [equilibrate](@ref Reaktoro::equilibrate) function use before to perform the equilibrium calculation.

--------------------------------

Using %Reaktoro with GEMS and PHREEQC
------------------------------------

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
<a href="https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/" target="_blank"> <img src="img/phreeqc-logo.jpg" width="100%"/></a>
</div>
<div class="col-md-1"></div>
</div>
</div>
@endhtmlonly

Because many users rely on specific thermodynamic models (e.g., models for activity coefficients, fugacity coefficients) and databases (e.g., PSI/NAGRA, PHREEQC databases) supported in these two software, two classes have been implemented in Reaktoro to permit the combined use of those models and databases with Reaktoro's numerical methods. These classes are [Gems](@ref Reaktoro::Gems) and [Phreeqc](@ref Reaktoro::Phreeqc).

The code below briefly (and incompletely) demonstrate how GEMS thermodynamic models and databases can be used with Reaktoro. It assumes a file named `project.lst` was output from the [GEM-Selektor][GEMS] application, which is then used to create an object of class [Gems](@ref Reaktoro::Gems). This object is then used to initialize the most important class in Reaktoro: the [ChemicalSystem](@ref Reaktoro::ChemicalSystem) class. Once this is done, chemical equilibrium and kinetics problems can be defined using this [ChemicalSystem](@ref Reaktoro::ChemicalSystem) object, like shown in the previous examples.

~~~{.cpp}
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Use exported file in GEM-Selektor to initialize an object of class Gems.
    Gems gems;("project.lst");

    // Use Gems object *gems* to construct a ChemicalSystem object, *system*,
    // that encapsulates information of the phases, species, and elements that
    // compose the chemical system, and associated thermodynamic models
    // specified in GEM-Selektor.
    ChemicalSystem system = gems;

    // Use Gems object to construct a ChemicalState object that contains
    // the temperature, pressure, and amounts of species produced by GEMS.
    ChemicalState state = gems.state(system);

    // Define here your chemical equilibrium and kinetics problems using
    // the just created ChemicalSystem object `system`.
    // These will be solved using Reaktoro's numerical methods.
    // GEMS will be used to calculate thermodynamic properties
    // of species (e.g., activities, standard chemical potentials) and
    // phases (e.g., molar volume, heat capacity, etc.).
    // The example below assumes the temperature used in GEMS was 25 °C,
    // and we now use Reaktoro's chemical equilibrium method to calculate
    // the new equilibrium state when temperature is changed to 50 °C.
    state.setTemperature(50, "celsius");

    // Perform the equilibrium calculation and update *state*.
    equilibrate(state);

    // Output the new chemical state to a file.
    state.output("state.txt");
}
~~~

~~~{.cpp}
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Phreeqc phreeqc;
    phreeqc.load("some-phreeqc-database-file.dat");
    phreeqc.execute("some-phreeqc-input-script.dat");

    // Use Phreeqc object to construct a ChemicalSystem object
    ChemicalSystem system = phreeqc;

    // From now on, use the `system`, an object of ChemicalSystem, to
    // perform chemical equilibrium and kinetics calculations using Reaktoro's numerical methods!
}
~~~

--------------------------------

Getting Started
---------------

@subpage installation Installation

@subpage tutorials Tutorials

@subpage about About



[PHREEQC]: http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
[GEMS]: http://gems.web.psi.ch/