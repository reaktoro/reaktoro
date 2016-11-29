# Tutorials {#tutorials}

We briefly show here some basic usage of Reaktoro for performing multiphase chemical equilibrium and kinetic calculations. More advanced use and customization is present in its user manual (work in progress!).

### Chemical equilibrium calculations
Reaktoro contains methods for chemical equilibrium calculations based on either the *Gibbs energy minimization* (GEM) approach or the *law of mass-action* (LMA) approach. In this section we describe step-by-step its use for performing chemical equilibrium calculations. 

#### A basic equilibrium calculation
In the code snippet below we show how the C++ interface of Reaktoro can be used to:

1. initialize a thermodynamic database;
2. specify the phases, and their species, composing the chemical system;
3. specify the equilibrium conditions for the calculation; and
4. perform the equilibrium calculation using a default method.

~~~{.cpp}
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Database db("supcrt98.xml");

    ChemicalEditor editor(db);
    editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- HCO3- CO2(aq) CO3--");
    editor.addGaseousPhase("H2O(g) CO2(g)");
    editor.addMineralPhase("Halite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.setTemperature(60, "celsius");
    problem.setPressure(300, "bar");
    problem.add("H2O", 1, "kg");
    problem.add("CO2", 100, "g");
    problem.add("NaCl", 1, "mol");

    ChemicalState state = equilibrate(problem);

    std::cout << state << std::endl;
}
~~~

The first two lines:

~~~{.cpp}
#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;
~~~

include the main Reaktoro header file: `Reaktoro.hpp`. By doing this, the application has access to all its classes and methods. The second line above is needed for convenience reasons: it eliminates the need to explicitly specify the namespace of Reaktoro components. Without it, we would need to write `Reaktoro::Database`, `Reaktoro::ChemicalSystem`,  `Reaktoro::EquilibriumProblem`, and so forth, which is a lot more verbose.

The equilibrium calculation uses the SUPCRT database together with the revised  *Helgeson-Kirkham-Flowers* (HKF) equations of state for the calculation of standard thermodynamic properties of aqueous, gaseous, and mineral species at temperatures 0 to 1000 °C and pressures 1 to 5000 bar. 

~~~{.cpp}
Database db("supcrt98.xml");
~~~

The line above is needed to initialize a `Database` object using a database file `supcrt98.xml` that should be found **at the same directory from where the application is executed!** Luckily, Reaktoro maintains a few built-in databases, including `supcrt98.xml`, whose files need not be found along with the application. Thus, if no file `supcrt98.xml` is found, the `Database` object will be initialized with the built-in database. 

The table below shows the current available built-in databases in Reaktoro and their description.

| Built-in Database | Description 
|-|-
| `supcrt98.xml` | Derived from the [SUPCRT][supcrt] database file `slop98.dat`, **without** organic species.
| `supcrt07.xml` | Derived from the [SUPCRT][supcrt] database file `slop07.dat`, **without** organic species.
| `supcrt98-organics.xml` | Derived from the [SUPCRT][supcrt] database file `slop98.dat`, **with** organic species.
| `supcrt07-organics.xml` | Derived from the [SUPCRT][supcrt] database file `slop07.dat`, **with** organic species.

Once the `Database` object has been initialized, one can use it to define the chemical system. For this, it is convenient to use the `ChemicalEditor` class, which currently permits the specification of aqueous, gaseous, and mineral phases, as well as specifying temperature and pressure interpolation points for the standard thermodynamic properties, and configuring the equations of state (e.g., HKF, Pitzer, Peng-Robinson, and many others) for calculation of activity/fugacity coefficients of the species. In the lines below we use the `ChemicalEditor` class to define a chemical system composed by an aqueous phase, a gaseous phase, and two pure-mineral phases. 

~~~{.cpp}
ChemicalEditor editor(db);
editor.addAqueousPhase("H2O(l) H+ OH- Na+ Cl- Ca++ HCO3- CO2(aq) CO3--");
editor.addGaseousPhase("H2O(g) CO2(g)");
editor.addMineralPhase("Halite");
editor.addMineralPhase("Calcite");
~~~
The names of the species listed above, **separated by spaces**, must be found in the provided database file `supcrt98.xml`, otherwise an exception will be thrown. Note that there can only be one aqueous and one gaseous phase in the chemical system. Calling methods `addAqueousPhase` and `addGaseousPhase` more than once will simply replace the definition of those phases. However, method `addMineralPhase` can be called as many times as there are mineral phases in the system. If the mineral phase is a solid solution, then specify the mineral end-members like it was done in `addAqueousPhase` and `addGaseousPhase`. Note that a chemical system does not necessarily need to have an aqueous phase, or a gaseous phase, or mineral phases. Choose the combination of phases that describes your problem.

After the chemical system has been defined using class `ChemicalEditor`, it is now time to initialize an instance of `ChemicalSystem` class:

~~~{.cpp}
ChemicalSystem system(editor);
~~~

The `ChemicalSystem` class is one of the most important classes in Reaktoro. It is the class used to computationally represent a chemical system, with all its phases, species, and elements. It is also the class used for computation of thermodynamic properties of phases and species, such as activities, chemical potentials, standard Gibbs energies, enthalpies, phase molar volumes, densities, and many others. Many classes in Reaktoro require an instance of `ChemicalSystem` for their initialization, since any chemical calculation needs to know the composition of the chemical system and the equations of state describing the non-ideal behavior of the phases.

Reaktoro provides the class `EquilibriumProblem` for convenient description of equilibrium conditions. Using this class allows one to set the temperature and pressure at equilibrium, and a recipe that describes a mixture of substances and their amounts, which can be seen as initial conditions for the equilibrium calculation:

~~~{.cpp}
EquilibriumProblem problem(system);
problem.setTemperature(60, "celsius");
problem.setPressure(300, "bar");
problem.add("H2O", 1, "kg");
problem.add("CO2", 100, "g");
problem.add("NaCl", 1, "mol");
~~~
The units above can be changed, or even suppressed. If not provided, default units are used, such as K for temperatures, Pa for pressures, and mol for amounts. The `add` method in `EquilibriumProblem` supports both amount and mass units, such as `mmol`,  `umol`, `g`, `mg`, etc.

!!! note ""
    <i class="fa fa-info-circle"></i> **Note**: Make sure you use these units consistently. Using temperature units when setting pressure, for example, will throw an exception!

Once the equilibrium problem has been defined, it is now time to solve it. This can be done using the utility method `equilibrate`:

~~~{.cpp}
ChemicalState state = equilibrate(problem);
~~~
The line above uses the definition of the equilibrium problem stored in object `problem` to perform the equilibrium calculation. The result of the calculation is the object `state`, an instance of  `ChemicalState` class, which is used to store the chemical state (i.e., the temperature, pressure, and molar amounts of all species) of the system at prescribed equilibrium conditions. The `ChemicalState` class contains also methods for querying thermodynamic properties of the system.

Finally, we output the chemical state of the system to the standard output using:

~~~{.cpp}
std::cout << state << std::endl;
~~~

This will output tables describing the chemical state of the system. For example, the molar amounts, molar fractions, activities, activity coefficients, and chemical potentials of the species. The molar volumes of the phases, the amounts of each element in the phases, and also the total phase molar amounts. The result of the above equilibrium problem can be seen in this [figure](../img/demo-equilibrium1-table.png).

#### Reaction path calculation: equilibrium-controlled reaction path

Consider two different chemical states in equilibrium: an *initial state* and a *final state*. These states can have different temperatures, pressures, and/or molar amounts of elements. If we gradually adjust temperature, pressure, and elemental amounts in the system to bring the initial state to the final state, slowly enough so that **every intermediate state is in equilibrium**, the system would trace a path, which we call *reaction path*. 

Let's say we initially have 100 g of calcite (CaCO3) mixed with 1 kg of water. We want to see how the addition of hydrochloric acid (HCl), up to 1mmol, contributes to the dissolution of calcite. Thus, our initial and final states for a reaction path calculation can be described as follows:

| Initial state  | Final state    |
|----------------|----------------|
| 1 kg of H2O    | 1 kg of H2O    |
| 100 g of CaCO3 | 100 g of CaCO3 |
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

    ChemicalSystem system(editor);

    EquilibriumProblem problem1(system);
    problem1.setTemperature(30.0, "celsius");
    problem1.setPressure(1.0, "bar");
    problem1.add("H2O", 1, "kg");
    problem1.add("CaCO3", 100, "g");

    EquilibriumProblem problem2(system);
    problem2.setTemperature(30.0, "celsius");
    problem2.setPressure(1.0, "bar");
    problem2.add("H2O", 1, "kg");
    problem2.add("CaCO3", 100, "g");
    problem2.add("HCl", 1, "mmol");

    ChemicalState state1 = equilibrate(problem1);
    ChemicalState state2 = equilibrate(problem2);

    EquilibriumPath path(system);

    ChemicalPlot plot0 = path.plot();
    plot0.x("pH");
    plot0.y("molality element=Ca units=molal");
    plot0.xlabel("pH");
    plot0.ylabel("Concentration [molal]");
    plot0.legend("Ca");

    ChemicalPlot plot1 = path.plot();
    plot1.x("amount element=Cl units=mmol");
    plot1.y("pH");
    plot1.xlabel("HCl [mmol]");
    plot1.ylabel("pH");
    plot1.nolegend();

    ChemicalOutput output = path.output();
    output.header("HCl [mmol]; Ca [molal]; pH");
    output.data("amount element=Cl units=mmol; molality element=Ca; pH");
    output.file("result.txt");

    path.solve(state1, state2);
}
~~~
In the code above, two instances of class `EquilibriumProblem` are created: `problem1` describes the initial state, and `problem2` describes the final state. Two instances of class `ChemicalState` are then created to store the initial and final equilibrium states calculated by method `equilibrate`. 

Note that, differently from the previous code example, the object `editor` from class `ChemicalEditor` was not initialized with a given `Database` object. Instead, it was initialized using the default built-in database file `supcrt98.xml`. Also note that the aqueous species were not listed, but the chemical elements composing the phase. Using element names to define any phase results in that phase containing all species in the database that can be built from those elements. 

Once the initial and final equilibrium states have been calculated, it is now time to trace the reaction path between them, with each intermediate state in chemical equilibrium. For this, we use the class `EquilibriumPath` in Reaktoro. Note that its initialization requires a `ChemicalSystem` instance:

~~~{.cpp}
EquilibriumPath path(system);
~~~

The results of a reaction path calculation is better analyzed via plots. Before the method `EquilibriumPath::solve` is called, one can configure plots to be generated at real-time, during the calculation. These plots are generated by [Gnuplot](http://www.gnuplot.info/), so ensure it is installed in your system if you want to see these plots. There are two configured plots in the previous code:

~~~{.cpp}
EquilibriumPath path(system);

ChemicalPlot plot0 = path.plot();
plot0.x("pH");
plot0.y("molality element=Ca units=molal");
plot0.xlabel("pH");
plot0.ylabel("Concentration [molal]");
plot0.legend("Ca");

ChemicalPlot plot1 = path.plot();
plot1.x("amount element=Cl units=mmol");
plot1.y("pH");
plot1.xlabel("HCl [mmol]");
plot1.ylabel("pH");
plot1.nolegend();
~~~

The first plot, `plot0`, is configured to have the x-axis containing the pH of the solution and the y-axis the molality of element Ca. The second plot, `plot1`, has the amount of element Cl, in units of mmol, as the x-axis, while the y-axis is pH of the solution. Thus, this second graph shows how pH changes with addition of HCl, while the first graph shows how the concentration of Ca changes with pH. These plots are shown below:

<div style="text-align: center; padding-bottom: 15px;">
    <a href="../img/fig-equilibriumpath-calcite-hcl-molality-ca.png" 
        data-lightbox="group1" data-title="The molality of element Ca with increase in pH.">
            <img src="../img/fig-equilibriumpath-calcite-hcl-molality-ca.png" width="40%"></a>

    <a href="../img/fig-equilibriumpath-calcite-hcl-ph.png" 
        data-lightbox="group1" data-title="The pH of the solution with addition of HCl.">
            <img src="../img/fig-equilibriumpath-calcite-hcl-ph.png" width="40%"></a>
</div>

It is possible to plot more than one quantity in the same graph. For example, to plot the activity coefficient of species H+ and the mass of calcite in addition to the molality of element Ca in the first plot, one could could write instead:

\todo Find a way to hide `using namespace Reaktoro;` below

~~~{.cpp}
using namespace Reaktoro; 
ChemicalPlot plot2 = path.plot();
plot2.y("molality element=Ca units=molal;"
        "activityCoefficient species=H+;"
        "mass species=Calcite units=g");
~~~
\todo Fix this warning - add admonition style.
\warning <i class="fa fa-warning"></i> **Warning:** The quantities to be plotted must be separated by `;`. Space is used to separate the quantity (e.g., molality) from the entity that it acts upon (element, species, phase, or reaction) and from its units. The `=` sign should not be separated by spaces. Use `element=H` and not `element = H`.

If you want to output quantities during the calculation to a file or terminal, then use method `EquilibriumPath::output`, which returns an instance of class `ChemicalOutput`:

~~~{.cpp}
ChemicalOutput out = path.output();
out.header("Amount of Cl; Molality of Ca; pH");
out.data("amount element=Cl units=mmol; molality element=Ca; pH");
out.file("result.txt");
~~~

The method `ChemicalOutput::header` sets the header of the output table, with each column title separated by `;`. The method `ChemicalOutput::data` sets the quantities to be output, also separated by `;`. Finally, the method `ChemicalOutput::file` sets the name of the output file. To output the result directly to the standard output, use:

~~~{.cpp}
out.terminal(true);
~~~

Finally, after all plots and output have been configured, the equilibrium path can be calculated via:

~~~{.cpp}
path.solve(state1, state2);
~~~

### Chemical kinetics calculations

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
    Database db("databases/supcrt/supcrt98.xml");

    ChemicalEditor editor(db);

    editor.addAqueousPhase("H2O HCl CaCO3");
    editor.addMineralPhase("Calcite");

    editor.addMineralReaction("Calcite")
        .setEquation("Calcite = Ca++ + CO3--")
        .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
        .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
        .setSpecificSurfaceArea(10, "cm2/g");

    ChemicalSystem system(editor);
    ReactionSystem reactions(editor);

    EquilibriumProblem problem(system);
    problem.setPartition("inert = Calcite");
    problem.add("H2O", 1, "kg");
    problem.add("HCl", 1, "mmol");

    ChemicalState state0 = equilibrate(problem);

    state0.setSpeciesAmount("Calcite", 100, "g");

    KineticPath path(reactions);
    path.setPartition("kinetic = Calcite");

    ChemicalPlot plot1 = path.plot();
    plot1.x("t units=minute");
    plot1.y("molality element=Ca");
    plot1.xlabel("t [minute]");
    plot1.ylabel("Concentration [molal]");
    plot1.legend("Ca");

    ChemicalPlot plot2 = path.plot();
    plot2.x("t units=minute");
    plot2.y("amount species=Calcite units=g");
    plot2.xlabel("t [minute]");
    plot2.ylabel("Amount [g]");
    plot2.legend("Calcite");

    path.solve(state0, 0, 5, "minute");
}
~~~

[supcrt]: http://geopig.asu.edu/?q=tools