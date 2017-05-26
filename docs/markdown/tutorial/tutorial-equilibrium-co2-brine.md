# Chemical equilibrium calculation for a H@sub{2}O–NaCl–CO@sub{2} system {#tutorial-equilibrium-co2-brine}

The code below calculates the equilibrium state for a H@sub{2}O–NaCl–CO@sub{2} system when 1 kg of H@sub{2}O, 100 g of CO@sub{2}, and 1 mol of NaCl are mixed at temperature 60 °C and pressure 300 bar:

~~~{.cpp}
    #include <Reaktoro/Reaktoro.hpp>
    using namespace Reaktoro;

    int main()
    {
        Database db;("supcrt98.xml");

        ChemicalEditor editor;(db);
        editor.addAqueousPhase("H2O NaCl CO2");
        editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
        editor.addMineralPhase("Halite");

        ChemicalSystem system;(editor);

        EquilibriumProblem problem;(system);
        problem.setTemperature(60, "celsius");
        problem.setPressure(300, "bar");
        problem.add("H2O", 1, "kg");
        problem.add("CO2", 100, "g");
        problem.add("NaCl", 0.1, "mol");

        ChemicalState state = equilibrate(problem);

        state.output("state.txt");
    }
~~~

@htmlonly
<a href="img/result-demo-equilibrium-co2-brine.txt" target="_blank"
    <button type="button" class="btn btn-primary">Check the Output »</button>
</a>
@endhtmlonly

The output text file contains details about the equilibrium state of the defined chemical system for the given equilibrium conditions. The output will include, for example, the *amounts*, *masses*, *mole fractions*, *activities*, *activity coefficients*, and *chemical potentials* of the chemical species in each phase. It will also list properties of each phase, such as *density*, *molar volume*, *volume fraction*, as well specific properties of some phases (e.g., *ionic strengh*, *pH*, *pe* for the aqueous phase).

The first line:

~~~{.cpp}
    #include <Reaktoro/Reaktoro.hpp>
~~~

include the main %Reaktoro header file, `Reaktoro.hpp`. By doing this, the application has access to all %Reaktoro's classes and methods. 

The second line:

~~~{.cpp}
    using namespace Reaktoro;
~~~

is needed for convenience reasons. It eliminates the need to explicitly specify the %Reaktoro namespace for each of its components. Without it, we would need to write Reaktoro::Database, Reaktoro::ChemicalSystem, Reaktoro::EquilibriumProblem, and so forth.

The next line of code, inside the `main` function:

~~~{.cpp}
using namespace Reaktoro; {delete}
Database db;("supcrt98.xml");
~~~

initializes a [Database](@ref Reaktoro::Database) object using the thermodynamic database file [supcrt98.xml]. This file was produced by converting the SUPCRT92@sup{@cite Johnson1992} database file [slop98.dat] into XML format. Read @ref thermodynamic-databases for more details about this and other databases supported by %Reaktoro.

Once the [Database](@ref Reaktoro::Database) object has been initialized, one can use it to construct an object of class [ChemicalSystem](@ref Reaktoro::ChemicalSystem). This is one of the most important classes in %Reaktoro, and it is the class that represents a chemical system as a collection of phases, with each phase containing one or more chemical species. 

The lines of code:

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalEditor editor;(db);
editor.addAqueousPhase("H2O NaCl CO2");
editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
editor.addMineralPhase("Halite");
~~~

use class [ChemicalEditor](@ref Reaktoro::ChemicalEditor) to define a chemical system with three phases: an *aqueous phase*, a *gaseous phase*, and a *mineral phase*. 

Defining a phase with more than one species requires either a *string containing multiple substance names*, such as in the call to [addAqueousPhase](@ref Reaktoro::ChemicalEditor::addAqueousPhase), or a *vector with species names*, such as in the call to [addGaseousPhase](@ref Reaktoro::ChemicalEditor::addGaseousPhase). When a string containing multiple substances are given, e.g., `"H2O NaCl CO2"`, the chemical species are selected automatically from the database based on the list of chemical elements that compose those substances. For example, the aqueous phase above will contain all aqueous species in the specified database that can be formed from elements H, O, Na, Cl, and C. Thus, one could alternatively have defined the previous aqueous phase using:

~~~{.cpp}
editor.addAqueousPhase("H O Na Cl C");
~~~

When a vector with species names are given instead, then the phase is constructed with only those species. Because of this, the names of the species must be found in the specified database file, `supcrt98.xml`, otherwise an exception will be thrown. Finally, defining a phase with only one chemical species can be done by either providing the species name in a string, as in the call to [addMineralPhase](@ref Reaktoro::ChemicalEditor::addMineralPhase) above, or using a vector of species names with only one name:

~~~{.cpp}
editor.addMineralPhase({"Halite"});
~~~

After the chemical system has been defined using class [ChemicalEditor](@ref Reaktoro::ChemicalEditor), it is now time to create an object of class [ChemicalSystem](@ref Reaktoro::ChemicalSystem):

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalSystem system(editor);
~~~

The [ChemicalSystem](@ref Reaktoro::ChemicalSystem) class is one of the most important classes in Reaktoro. It is the class used to computationally represent a chemical system, with all its phases, species, and elements. It is also the class used for computation of thermodynamic properties of phases and species, such as activities, chemical potentials, standard Gibbs energies, enthalpies, phase molar volumes, densities, and many others. Many classes in %Reaktoro require an instance of [ChemicalSystem](@ref Reaktoro::ChemicalSystem) for their initialization, since any chemical calculation needs to know the definition of the chemical system and the thermodynamic models describing the non-ideal behavior of the phases. Read @ref defining-chemical-systems for more information on the various ways one can define the chemical system and how thermodynamic models can be specified for each phase.

%Reaktoro provides the class [EquilibriumProblem](@ref Reaktoro::EquilibriumProblem) for convenient description of equilibrium conditions. Using this class allows one to set the temperature and pressure at equilibrium, and a recipe that describes a mixture of substances and their amounts, which can be seen as initial conditions for the equilibrium calculation:

~~~{.cpp}
using namespace Reaktoro; {delete}
EquilibriumProblem problem;(system);
problem.setTemperature(60, "celsius");
problem.setPressure(300, "bar");
problem.add("H2O", 1, "kg");
problem.add("CO2", 100, "g");
problem.add("NaCl", 0.1, "mol");
~~~

The units above can be changed, or even suppressed. If not provided, default units are used, such as K for temperatures, Pa for pressures, and mol for amounts. The `add` method in [EquilibriumProblem](@ref Reaktoro::EquilibriumProblem) supports both amount and mass units, such as `mmol`,  `umol`, `g`, `mg`, etc.

@warning Make sure units are used consistently. Using temperature units when setting pressure, for example, will result in an error.

Once the equilibrium problem has been defined, it is now time to solve it. This can be done using the utility method [equilibrate](@ref Reaktoro::equilibrate):

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalState state = equilibrate(problem);
~~~

The line above uses the definition of the equilibrium problem stored in the object `problem` to perform the equilibrium calculation. The result of the calculation is the object `state`, an instance of  [ChemicalState](@ref Reaktoro::ChemicalState) class, which is used to store the chemical state (i.e., the temperature, pressure, and molar amounts of all species) of the system at prescribed equilibrium conditions. The [ChemicalState](@ref Reaktoro::ChemicalState) class also provides methods for querying thermodynamic properties of the system (see method [ChemicalState::properties](@ref Reaktoro::ChemicalState::properties)).

@note Method [equilibrate](@ref Reaktoro::equilibrate) is not the optimal method for performing a sequence of equilibrium calculations (e.g., when coupling Reaktoro with other codes for simulating fluid flow, species transport, etc.). In situations where many equilibrium calculations need to be performed and good initial guesses are available each time (e.g., the equilibrium state at a previous time step serving as an initial guess for the new equilibrium calculation), use class [EquilibriumSolver](@ref Reaktoro::EquilibriumSolver).

Finally, we output the chemical state of the system to a file named `state.txt` using:

~~~{.cpp}
state.output("state.txt");
~~~

[supcrt]: http://geopig.asu.edu/?q=tools
[supcrt98.xml]: databases/supcrt/supcrt98.xml
[supcrt07.xml]: databases/supcrt/supcrt07.xml
[slop98.dat]: databases/supcrt/slop98.dat
[slop07.dat]: databases/supcrt/slop07.dat

