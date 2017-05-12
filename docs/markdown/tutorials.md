# Tutorials {#tutorials}

We briefly show here some basic usage of Reaktoro for performing multiphase chemical equilibrium and kinetic calculations.

## Chemical equilibrium calculations
Reaktoro contains methods for chemical equilibrium calculations based on either the *Gibbs energy minimization* (GEM) approach or the *law of mass-action* (LMA) approach. In this section, we describe step-by-step how it can be used for performing chemical equilibrium calculations.

### A basic equilibrium calculation
In the code snippet below we show how the C++ interface of Reaktoro can be used to:

1. initialize a thermodynamic database;
2. specify the phases, and their species, that compose the chemical system of interest;
3. specify the equilibrium conditions for the calculation; and
4. perform the equilibrium calculation using a default Gibbs energy minimization method.

~~~{.cpp}
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database db;("supcrt98.xml");

    ChemicalEditor editor;(db);
    editor.addAqueousPhase({"H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"});
    editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
    editor.addMineralPhase("Halite");

    ChemicalSystem system;(editor);

    EquilibriumProblem problem;(system);
    problem.setTemperature(60, "celsius");
    problem.setPressure(300, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CO2", 100, "g");
    problem.add("NaCl", 1, "mol");

    EquilibriumState state = equilibrate(problem);

    std::cout << state << std::endl;
}
~~~

The first two lines:

~~~{.cpp}
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;
~~~

include the main Reaktoro header file: Reaktoro.hpp. By doing this, the application has access to all its classes and methods. The second line above is needed for convenience reasons: it eliminates the need to explicitly specify the namespace of Reaktoro components. Without it, we would need to write Reaktoro::Database, Reaktoro::ChemicalSystem,  Reaktoro::EquilibriumProblem, and so forth, which is a lot more verbose.

The equilibrium calculation uses the SUPCRT database together with the revised  *Helgeson-Kirkham-Flowers* (HKF) equations of state for the calculation of standard thermodynamic properties of aqueous, gaseous, and mineral species at temperatures 0 to 1000 °C and pressures 1 to 5000 bar.

~~~{.cpp}
using namespace Reaktoro; {delete}
Database db;("supcrt98.xml");
~~~

The line above is needed to initialize a [Database](@ref Reaktoro::Database) object using a database file `supcrt98.xml` that should be found **at the same directory from where the application is executed!** For convenience, Reaktoro also maintains a few built-in databases, including `supcrt98.xml`, whose files need not be found along with the application. Thus, if no file `supcrt98.xml` is found, the [Database](@ref Reaktoro::Database) object will be initialized with the built-in database instead.

The table below shows the current available built-in databases in Reaktoro and their description.

| Built-in Database | Description
|-|-
| `supcrt98.xml` | Derived from the [SUPCRT][supcrt] database file `slop98.dat`, **without** organic species.
| `supcrt07.xml` | Derived from the [SUPCRT][supcrt] database file `slop07.dat`, **without** organic species.
| `supcrt98-organics.xml` | Derived from the [SUPCRT][supcrt] database file `slop98.dat`, **with** organic species.
| `supcrt07-organics.xml` | Derived from the [SUPCRT][supcrt] database file `slop07.dat`, **with** organic species.

Once the [Database](@ref Reaktoro::Database) object has been initialized, one can use it to define the chemical system. For this, it is convenient to use the [ChemicalEditor](@ref Reaktoro::ChemicalEditor) class, which currently permits the specification of aqueous, gaseous, and mineral phases, as well as specifying temperature and pressure interpolation points for the standard thermodynamic properties, and configuring the equations of state (e.g., HKF, Pitzer, Peng-Robinson, and many others) for calculation of activity/fugacity coefficients of the species. In the lines below, we use the [ChemicalEditor](@ref Reaktoro::ChemicalEditor) class to define a chemical system composed by an aqueous phase, a gaseous phase, and a pure mineral phase.

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalEditor editor;(db);
editor.addAqueousPhase({"H2O(l)", "H+", "OH-", "Na+", "Cl-", "HCO3-", "CO2(aq)", "CO3--"});
editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
editor.addMineralPhase("Halite");
~~~

Defining a phase with more than one species requires a list of species names, like seen for the [addAqueousPhase](@ref Reaktoro::ChemicalEditor::addAqueousPhase) and [addGaseousPhase](@ref Reaktoro::ChemicalEditor::addGaseousPhase) method calls, while for a phase with a single species, like the pure mineral phase halite, with chemical formula NaCl(s), a string with a single name suffices.

@warning The names of the species listed above must be found in the provided database file `supcrt98.xml`, otherwise an exception will be thrown. We will later see an alternative approach in which the species in a phase can be **automatically** selected from the database from a given list of chemical elements or substance names.

@note Currently, there can only be one aqueous and one gaseous phase in the chemical system. Calling methods [addAqueousPhase](@ref Reaktoro::ChemicalEditor::addAqueousPhase) and [addGaseousPhase](@ref Reaktoro::ChemicalEditor::addGaseousPhase) more than once will simply replace the definition of those phases. However, method [addMineralPhase](@ref Reaktoro::ChemicalEditor::addMineralPhase) can be called as many times as there are mineral phases in the system. If the mineral phase is a solid solution, then specify the mineral end-members like it was done in [addAqueousPhase](@ref Reaktoro::ChemicalEditor::addAqueousPhase) and [addGaseousPhase](@ref Reaktoro::ChemicalEditor::addGaseousPhase). Note that a chemical system does not necessarily need to have an aqueous phase, or a gaseous phase, or mineral phases. Choose the combination of phases that describes your problem.

After the chemical system has been defined using class [ChemicalEditor](@ref Reaktoro::ChemicalEditor), it is now time to initialize an instance of [ChemicalSystem](@ref Reaktoro::ChemicalSystem) class:

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalSystem system(editor);
~~~

The [ChemicalSystem](@ref Reaktoro::ChemicalSystem) class is one of the most important classes in Reaktoro. It is the class used to computationally represent a chemical system, with all its phases, species, and elements. It is also the class used for computation of thermodynamic properties of phases and species, such as activities, chemical potentials, standard Gibbs energies, enthalpies, phase molar volumes, densities, and many others. Many classes in Reaktoro require an instance of [ChemicalSystem](@ref Reaktoro::ChemicalSystem) for their initialization, since any chemical calculation needs to know the composition of the chemical system and the thermodynamic models describing the non-ideal behavior of the phases.

Reaktoro provides the class [EquilibriumProblem](@ref Reaktoro::EquilibriumProblem) for convenient description of equilibrium conditions. Using this class allows one to set the temperature and pressure at equilibrium, and a recipe that describes a mixture of substances and their amounts, which can be seen as initial conditions for the equilibrium calculation:

~~~{.cpp}
using namespace Reaktoro; {delete}
EquilibriumProblem problem;(system);
problem.setTemperature(60, "celsius");
problem.setPressure(300, "bar");
problem.add("H2O", 1, "kg");
problem.add("CO2", 100, "g");
problem.add("NaCl", 1, "mol");
~~~

The units above can be changed, or even suppressed. If not provided, default units are used, such as K for temperatures, Pa for pressures, and mol for amounts. The `add` method in [EquilibriumProblem](@ref Reaktoro::EquilibriumProblem) supports both amount and mass units, such as `mmol`,  `umol`, `g`, `mg`, etc.

@warning Make sure units are used consistently. Using temperature units when setting pressure, for example, will result in an error.

Once the equilibrium problem has been defined, it is now time to solve it. This can be done using the utility method [equilibrate](@ref Reaktoro::equilibrate):

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalState state = equilibrate(problem);
~~~

The line above uses the definition of the equilibrium problem stored in the object `problem` to perform the equilibrium calculation. The result of the calculation is the object `state`, an instance of  [EquilibriumState](@ref Reaktoro::EquilibriumState) class, which is used to store the chemical state (i.e., the temperature, pressure, and molar amounts of all species) of the system at prescribed equilibrium conditions. The [EquilibriumState](@ref Reaktoro::EquilibriumState) class also provides methods for querying thermodynamic properties of the system (see method [EquilibriumState::properties](@ref Reaktoro::EquilibriumState::properties).

@note Method [equilibrate](@ref Reaktoro::equilibrate) is not the optimal method for performing a sequence of equilibrium calculations (e.g., when coupling Reaktoro with other codes for simulating fluid flow, species transport, etc.). In situations where many equilibrium calculations need to be performed and good initial guesses are available each time (e.g., the equilibrium state at a previous time step serving as an initial guess for the new equilibrium calculation), use class [EquilibriumSolver](@ref Reaktoro::EquilibriumSolver).

Finally, we output the chemical state of the system to the standard output using:

~~~{.cpp}
std::cout << state << std::endl;
~~~

This will output tables describing the chemical state of the system. For example, the molar amounts, molar fractions, activities, activity coefficients, and chemical potentials of the species. The molar volumes of the phases, the amounts of each element in the phases, and also the total phase molar amounts. The result of the above equilibrium problem can be seen in this [figure](../img/demo-equilibrium1-table.png).

### Reaction path calculation: equilibrium-controlled reaction path

Consider two different chemical states in equilibrium: an *initial state* and a *final state*. These states can have different temperatures, pressures, and/or molar amounts of elements. If we gradually adjust temperature, pressure, and elemental amounts in the system to bring the initial state to the final state, slowly enough so that **every intermediate state is in equilibrium**, the system would trace a path, which we call *reaction path*.

Let's say we initially have 1 g of calcite (CaCO3) mixed with 1 kg of water. We want to see how the addition of hydrochloric acid (HCl), up to 1 mmol, contributes to the dissolution of calcite. Thus, our initial and final states for a reaction path calculation can be described as follows:

| Initial state  | Final state    |
|----------------|----------------|
| 1 kg of H2O    | 1 kg of H2O    |
| 1 g of CaCO3   | 1 g of CaCO3   |
|                | 1 mmol of HCl  |

We show below the code for doing this reaction path calculation using Reaktoro, followed by comments of each newly introduced C++ component:

~~~{.cpp}
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    ChemicalEditor editor;
    editor.addAqueousPhase("H O Ca C Cl");
    editor.addMineralPhase("Calcite");

    ChemicalSystem system;(editor);

    EquilibriumProblem problem1;(system);
    problem1.setTemperature(30.0, "celsius");
    problem1.setPressure(1.0, "bar");
    problem1.add("H2O", 1, "kg");
    problem1.add("CaCO3", 1, "g");

    EquilibriumProblem problem2;(system);
    problem2.setTemperature(30.0, "celsius");
    problem2.setPressure(1.0, "bar");
    problem2.add("H2O", 1, "kg");
    problem2.add("CaCO3", 1, "g");
    problem2.add("HCl", 1, "mmol");

    EquilibriumState state1 = equilibrate(problem1);
    EquilibriumState state2 = equilibrate(problem2);

    EquilibriumPath path;(system);

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
~~~

In the code above, two instances of class [EquilibriumProblem](@ref Reaktoro::EquilibriumProblem) are created: `problem1` describes the initial state, and `problem2` the final state. Two instances of class [ChemicalState](@ref Reaktoro::ChemicalState) are then created to store the initial and final equilibrium states calculated by method [equilibrate](@ref Reaktoro::equilibrate).

Note that, differently from the previous code example, the object `editor` from class [ChemicalEditor](@ref Reaktoro::ChemicalEditor) was not initialized with a given [Database](@ref Reaktoro::Database) object. Instead, it was initialized using the default built-in database file `supcrt98.xml`. Also note that the aqueous species were not listed, but the chemical elements composing the phase. When chemical element names are specified during the creation of a phase, like in:

~~~{.cpp}
using namespace Reaktoro; {delete}
editor.addAqueousPhase("H O Ca C Cl");
~~~

class [ChemicalEditor](@ref Reaktoro::ChemicalEditor) searches for all species in the database that can be formed by those elements. Only species corresponding to the phase type being created is selected (e.g., only aqueous species are searched in the above case).

Once the initial and final equilibrium states have been calculated, it is now time to trace the reaction path between them, with each intermediate state in chemical equilibrium. For this, we use the class [EquilibriumPath](@ref Reaktoro::EquilibriumPath). Note that its initialization requires a [ChemicalSystem](@ref Reaktoro::ChemicalSystem) instance:

~~~{.cpp}
using namespace Reaktoro; {delete}
EquilibriumPath path(system);
~~~

A reaction path calculation is better analyzed using plots. Before the method [EquilibriumPath::solve](@ref Reaktoro::EquilibriumPath::solve) is called, one can configure plots to be generated during the calculation. These plots are generated by [Gnuplot](http://www.gnuplot.info/), so ensure it is installed in your system to be able to see these plots. There are four configured plots in the above equilibrium path calculation. The first plot is defined as:

~~~{.cpp}
using namespace Reaktoro; {delete}
EquilibriumPath path;(system); {delete}
ChemicalPlot plot1 = path.plot();
plot1.x("elementAmount(Cl units=mmol)");
plot1.y("pH");
plot1.xlabel("HCl [mmol]");
plot1.ylabel("pH");
plot1.showlegend(false);
~~~

which sets the *x*-axis to the amount of element Cl, in units of mmol, and the *y*-axis to the pH of the aqueous phase, resulting in the following figure:

@htmlonly
<div class="image">
<a href="img/equilibrium-path-calcite-hcl-1.svg"
    data-fancybox="gallery1"
    data-caption="Tracing the change in pH with the gradual addition of HCl.">
    <img src="img/equilibrium-path-calcite-hcl-1.svg" width=400></a>
</div>
@endhtmlonly

The second plot is defined as:

~~~{.cpp}
using namespace Reaktoro; {delete}
EquilibriumPath path;(system); {delete}
ChemicalPlot plot2 = path.plot();
plot2.x("elementAmount(Cl units=mmol)");
plot2.y("Ca", "elementMolality(Ca units=mmolal)");
plot2.xlabel("HCl [mmol]");
plot2.ylabel("Concentration [mmolal]");
plot2.legend("right center");
~~~

which now sets the *x*-axis to the pH of the aqueous phase and the *y*-axis to the molality of element Ca, i.e., the molar amount of Ca **in the aqueous phase**, divided by the mass of solvent water H2O(l). This results in the following figure:

@htmlonly
<div class="image">
<a href="img/equilibrium-path-calcite-hcl-2.svg"
    data-fancybox="gallery1"
    data-caption="Tracing the change in concentration of element Ca in the aqueous solution with the gradual addition of HCl.">
    <img src="img/equilibrium-path-calcite-hcl-2.svg" width=400></a>
</div>
@endhtmlonly

The third figure is produced with the following code:

~~~{.cpp}
using namespace Reaktoro; {delete}
EquilibriumPath path;(system); {delete}
ChemicalPlot plot3 = path.plot();
plot3.x("elementAmount(Cl units=mmol)");
plot3.y("CO2(aq)", "speciesMolality(CO2(aq) units=mmolal)");
plot3.y("CO3--", "speciesMolality(CO3-- units=mmolal)");
plot3.xlabel("HCl [mmol]");
plot3.ylabel("Concentration [mmolal]");
plot3.legend("right top");
~~~

which also sets the *x*-axis to pH, but the *y*-axis now contains two plotted quantities: the molality of species CO2(aq) and the molality of species CO3--, both in units of mmolal (i.e., mmol/kgH2O). This produces the following figure:

@htmlonly
<div class="image">
<a href="img/equilibrium-path-calcite-hcl-3.svg"
    data-fancybox="gallery1"
    data-caption="Tracing the concentrations of species CO2(aq) and CO3-- with the gradual addition of HCl.">
    <img src="img/equilibrium-path-calcite-hcl-3.svg" width=400></a>
</div>
@endhtmlonly

The fourth and last figure finally plots how the mass of calcite (or calcium carbonate) changes with the addition of HCl in the system:

~~~{.cpp}
using namespace Reaktoro; {delete}
EquilibriumPath path;(system); {delete}
ChemicalPlot plot4 = path.plot();
plot4.x("elementAmount(Cl units=mmol)");
plot4.y("Calcite", "speciesMass(Calcite units=g)");
plot4.xlabel("HCl [mmol]");
plot4.ylabel("Mass [g]");
~~~

producing the following figure:

@htmlonly
<div class="image">
<a href="img/equilibrium-path-calcite-hcl-4.svg"
    data-fancybox="gallery1"
    data-caption="Tracing the dissolved mass of calcite with the gradual addition of HCl.">
    <img src="img/equilibrium-path-calcite-hcl-4.svg" width=400></a>
</div>
@endhtmlonly

@note Check class [ChemicalQuantity](@ref Reaktoro::ChemicalQuantity) for a list of supported quantity names, their default units, and how they can be used in both [ChemicalPlot](@ref Reaktoro::ChemicalPlot) and [ChemicalOutput](@ref Reaktoro::ChemicalOutput) classes.

@warning The equal sign (=) used to specify the units (e.g., `units=mmolal`, `units=celsius`) should not be separated by spaces.

To output quantities to a file or terminal during the calculation, use method [EquilibriumPath::output](@ref Reaktoro::EquilibriumPath::output), which returns an instance of class [ChemicalOutput](@ref Reaktoro::ChemicalOutput):

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalOutput output = path.output();
output.file("result.txt");
output.add("Cl [mmol]", "elementAmount(Cl units=mmol)");
output.add("Ca [mmolal]", "elementMolality(Ca units=mmolal)");
output.add("pH");
output.add("speciesMass(Calcite units=g)");
~~~

The method [ChemicalOutput::add](@ref Reaktoro::ChemicalOutput::add) adds a quantity to be output to the file `result.txt` (i.e., the file name specified in the call to method [ChemicalOutput::file](@ref Reaktoro::ChemicalOutput::file)). Each call to [ChemicalOutput::add](@ref Reaktoro::ChemicalOutput::add) results in a new column of data in the output file, like shown below:

~~~{.txt}
Cl [mmol]           Ca [mmolal]         pH                  speciesMass(Calcite units=g)
1.1e-16             0.134437            9.78558             0.986545            
6.3132e-06          0.134439            9.78557             0.986544            
9.55315e-05         0.134467            9.7853              0.986542            
0.000987715         0.134742            9.78268             0.986514            
0.00990955          0.137677            9.7559              0.98622             
0.029156            0.144655            9.69602             0.985522            
0.0484024           0.15261             9.63351             0.984726            
0.100331            0.179333            9.45629             0.982051            
0.152259            0.213487            9.27833             0.978633            
0.204187            0.253494            9.11218             0.974629            
0.256115            0.297405            8.96347             0.970234            
0.308043            0.343719            8.83241             0.965598            
0.397568            0.42662             8.6418              0.957301            
0.487093            0.511357            8.48573             0.94882             
0.531083            0.553286            8.41877             0.944623            
0.575074            0.595296            8.35705             0.940419            
0.605964            0.624812            8.31647             0.937465            
0.636854            0.654325            8.27795             0.934511            
0.68543             0.700704            8.22109             0.929869            
0.734007            0.747012            8.16829             0.925234            
0.797192            0.807094            8.10488             0.919221            
0.860377            0.866961            8.04663             0.913229            
0.922344            0.925434            7.99384             0.907376            
0.984311            0.983645            7.94481             0.90155             
1                   0.998339            7.93293             0.900079            
~~~

When two arguments are provided to method [ChemicalOutput::add](@ref Reaktoro::ChemicalOutput::add), the first one is a label used as the heading of the column of data in the output file, and the second argument is the name of the quantity to be output (e.g., `time`, `elementAmount(Cl)`, `ionicStrength`). When only one argument is provided, this single argument is both the label and the quantity name.

To output the result directly to the standard output, use:

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalOutput output; {delete}
output.terminal(true);
~~~

Finally, after all plots and output files have been configured, the equilibrium path can be calculated using:

~~~{.cpp}
using namespace Reaktoro; {delete}
EquilibriumPath path; {delete}
path.solve(state1, state2);
~~~

## 
We now finish by showing the Python code equivalent to the previous C++ code used for the equilibrium path calculation:

~~~{.cpp}
using namespace Reaktoro; {delete}
from reaktoro import *

editor = ChemicalEditor()
; {delete}
ChemicalEditor editor; {delete}
editor.addAqueousPhase("H O Ca C Cl")
; {delete}
editor.addMineralPhase("Calcite")
; {delete}

system = ChemicalSystem;(editor)

problem1 = EquilibriumProblem;(system)
problem1.setTemperature;(30.0, "celsius")
problem1.setPressure;(1.0, "bar")
problem1.add;("H2O", 1, "kg")
problem1.add;("CaCO3", 1, "g")

problem2 = EquilibriumProblem;(system)
problem2.setTemperature;(30.0, "celsius")
problem2.setPressure;(1.0, "bar")
problem2.add;("H2O", 1, "kg")
problem2.add;("CaCO3", 1, "g")
problem2.add;("HCl", 1, "mmol")

state1 = equilibrate(problem1)
state2 = equilibrate(problem2)

path = EquilibriumPath(system)
;{delete}

plot1 = path.plot()
plot1.x;("elementAmount(Cl units=mmol)")
plot1.y;("pH")
plot1.xlabel;("HCl [mmol]")
plot1.ylabel;("pH")
plot1.showlegend;(False)

plot2 = path.plot()
plot2.x("elementAmount(Cl units=mmol)")
plot2.y("Ca", "elementMolality(Ca units=mmolal)")
plot2.xlabel("HCl [mmol]")
plot2.ylabel("Concentration [mmolal]")
plot2.legend("right center")

plot3 = path.plot()
plot3.x("elementAmount(Cl units=mmol)")
plot3.y("CO2(aq)", "speciesMolality(CO2(aq) units=mmolal)")
plot3.y("CO3--", "speciesMolality(CO3-- units=mmolal)")
plot3.xlabel("HCl [mmol]")
plot3.ylabel("Concentration [mmolal]")
plot3.legend("right top")

plot4 = path.plot()
plot4.x("elementAmount(Cl units=mmol)")
plot4.y("Calcite", "speciesMass(Calcite units=g)")
plot4.xlabel("HCl [mmol]")
plot4.ylabel("Mass [g]")

output = path.output()
output.file("result.txt")
output.add("Cl [mmol]", "elementAmount(Cl units=mmol)")
output.add("Ca [mmolal]", "elementMolality(Ca units=mmolal)")
output.add("pH")
output.add("speciesMass(Calcite units=g)")

path.solve(state1, state2)
~~~

## Chemical kinetics calculations

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

    KineticState state0 = equilibrate(problem);

    state0.setSpeciesMass("Calcite", 100, "g");

    KineticPath path;(reactions);
    path.setPartition(partition);

    ChemicalPlot plot1 = path.plot();
    plot1.x("time(units=minute)");
    plot1.y("Ca", "elementMolality(Ca units=mmolal)");
    plot1.xlabel("Time [minute]");
    plot1.ylabel("Concentration [mmolal]");
    plot1.legend("right center");

    ChemicalPlot plot2 = path.plot();
    plot2.x("time(units=minute)");
    plot2.y("Calcite", "phaseMass(Calcite units=g)");
    plot2.xlabel("Time [minute]");
    plot2.ylabel("Mass [g]");

    path.solve(state0, 0, 5, "minute");
}

~~~

[supcrt]: http://geopig.asu.edu/?q=tools