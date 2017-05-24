# Defining the Chemical System {#defining-chemical-systems}

@tableofcontents

In Reaktoro, we need to define our chemical system in order to perform chemical reaction calculations.

# What is a chemical system? {#what-is-chemical-system}

A chemical system is a description of the phases of interest in the modeling problem and the chemical species that compose those phases. For example, when modeling the chemistry of an aqueous solution, only one phase should be enough, an *aqueous phase*. If one is interested in modeling the solubility of gases in an aqueous solution, then it makes sense to also define a *gaseous phase* with one or more gases. When modeling aqueous solutions and minerals, under a variety of temperature, pressure, and chemical conditions, it might make sense to define the chemical system with many *mineral phases*.

## Chemical and Thermodynamic Models {#chemical-thermo-models}

Specifying the phases and their species is not enough to fully describe a chemical system in the *computational sense*. Every phase in %Reaktoro has two associated models: a *thermodynamic model* and a *chemical model*. These denominations are not standard in the literature, but they are useful in the differentiation of two needed types of models for a phase.

A *thermodynamic model* is a model for the calculation of *standard thermodynamic properties* of the species. Examples include standard Gibbs energies, or standard chemical potentials, standard molar volumes, standard heat capacities, standard enthalpies, and so forth. These models are functions of *temperature* and *pressure* only. Currently, %Reaktoro natively supports only SUPCRT92 databases, which contains parameters for the revised  *Helgeson-Kirkham-Flowers* (HKF) equations of state for the calculation of standard thermodynamic properties for hundreds of aqueous species, at temperatures 0 to 1000 °C and pressures 1 to 5000 bar. The SUPCRT92 databases also contain Maier--Kelly coefficients for the calculation of standard thermodynamic properties of gases and minerals.

A *chemical model* is a model that describes the *non-ideal behavior* of phases. These models not only depend on temperature and pressure, like the thermodynamic models, but also on the amounts of the species in the phase. To be more precise, on the concentrations of these species, which can be calculated from the amounts of the species.

@todo Continue discussion on chemical and thermodynamic models.

## Deciding between chemical realism and computational efficiency {#realism-efficiency}

One can define a chemical system with many phases, each phase containing one or more species. This does not mean that all phases and their species exist at positive amounts! What it means is that the chemical calculations, equilibrium or kinetics, are capable of deciding if a phase in the chemical system should exist at positive amounts for some given conditions (e.g., temperature, pressure, overall composition).

By selecting as many phases as possible, with the possibilities constrained by the *thermodynamic database* being used, one can increase the confidence level of the estimated chemical states. Note, however, that accurate and realistic estimates depend on many more factors than just the selection of potential phases, such as  the *choice of thermodynamic models for non-ideal phases*. Furthermore, note that adding too many phases and species to the definition of the chemical system can result in *more computationally expensive* chemical calculations. In critical performance applications, such as when combining chemical reactions and fluid flow and species transport modeling, restricting the number of phases and species might be necessary for achieving feasible simulation times. The modeler is responsible to decide to which extent the number of phases and species can be compromised for efficiency reasons at the expense of chemical realism!


# Defining chemical systems in Reaktoro {#how-define}

Defining a chemical system in %Reaktoro can be done in different ways: all of them ultimately producing an object of class [ChemicalSystem]. In what follows, we demonstrate how such ChemicalSystem object can be created.

## Using class ChemicalEditor {#using-chemicaleditor}

%Reaktoro implements the class [ChemicalEditor] to facilitate the definition of a chemical system. To define a chemical system with only one aqueous phase, using the thermodynamic database [supcrt98.xml](@ref tutorial-database), one does:

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalEditor editor;("supcrt98.xml");
editor.addAqueousPhase({"H2O(l)", "H+", "OH-", "Na+", "Cl-"});

ChemicalSystem system = editor;
~~~

The aqueous phase is constructed with just five aqueous species. The species names must be written exactly as in the specified database, otherwise an error will occur! Once all phases, and their species have been specified in the `ChemicalEditor` object, one can now convert `editor` into an object of class `ChemicalSystem`, which can be done using the assignment operator `=`, which is programmed to work as a conversion operator:

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalSystem system = editor;
~~~

In Python, however, the equivalent code would be:

~~~{.cpp}
using namespace Reaktoro; {delete}
system = ChemicalSystem(editor)
~~~

In some occasions, it would be preferable to just specify a few compound or substance names, *not necessarily named as in the database*, and then let `ChemicalEditor` to automatically select all chemical species that could be formed out of the combination of those compounds. In this case, one no longer uses a list of species names, but a string containing all compound names that represents the phase. For the previous H2O--NaCl aqueous solution, the aqueous phase would be defined as follows:

~~~{.cpp}
using namespace Reaktoro; {delete}
editor.addAqueousPhase("H2O NaCl");
~~~

Neither `H2O` nor `NaCl` are species names in the `supcrt98.xml` database. Internally, a parser in %Reaktoro is used to decompose each compound or substance name into a list of chemical elements and their number of atoms. For example, parsing `H2O` would produce `{{"H", 2}, {"O", 1}}`. After this is done for all substance names, `ChemicalEditor` is then programmed to search in the database for all aqueous species that can be formed with these chemical elements. Because of this, one should expect many more chemical species in the aqueous phase than possibly necessary for the modeling application. If many chemical calculations are needed, so that computational time is critical in your application, then manually selecting the major chemical species might be a preferred option. This could be done by the modeler after a few analysis to detect which aqueous species are present in significant amounts in all or most possible conditions expected in the modeling.

To define a chemical system with a gaseous phase, say, for example, a H2O--NaCl--CO2 system, one does:

~~~{.cpp}
using namespace Reaktoro; {delete}
editor.addAqueousPhase("H2O NaCl CO2");
editor.addGaseousPhase("H2O CO2");
~~~

However, as mentioned above, this might result in many more species than needed, especially for the gaseous phase. For example, the resulting gaseous phase will most probably contain gases that might not be of interest for a H2O--NaCl--CO2 system, such as `CH4(g)`, `CO(g)`, `H2(g)`, `O2(g)`, etc. If water vapor, `H2O(g)`, and gaseous/supercritical carbon dioxide, `CO2(g)`, suffices, then the gaseous phase can be defined using:

~~~{.cpp}
using namespace Reaktoro; {delete}
editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
~~~

In most geochemical modeling applications, one or more mineral phases are needed. In most cases, these mineral phases are *pure mineral phases*, i.e., they contain only one mineral species. When they contain more than one, they are often called *solid solutions*. Defining a pure mineral phase or a solid solution phase is similar to defining any other phase type. The code below demonstrates the definition of a H2O--NaCl--CaCl2--MgCl2--CO2 system with pure mineral phases halite (NaCl), calcite (CaCO3), magnesite (MgCO3), and dolomite (CaMg(CO3)2):

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalEditor editor;("supcrt98.xml");
editor.addAqueousPhase("H2O NaCl CaCl2 MgCl2 CO2");
editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
editor.addMineralPhase("Halite");
editor.addMineralPhase("Calcite");
editor.addMineralPhase("Magnesite");
editor.addMineralPhase("Dolomite");

ChemicalSystem system = editor;
~~~

To define a solid solution phase, one would proceed similarly as done before for aqueous and gaseous phases:

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalEditor editor; {delete}
editor.addMineralPhase({"Calcite", "Magnesite"});
~~~

The code below demonstrate how an aqueous phase can be constructed and then have its default *chemical model* changed to the Debye--Huckel model:

~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalEditor editor; {delete}
AqueousPhase aqueousphase = editor.addAqueousPhase("H2O NaCl CO2");
aqueousphase.setChemicalModelDebyeHuckel();
~~~

In general, chemical models for aqueous phases Note that activity models are also needed for neutral species. If none are provided, then an ideal model is used in which their activity coefficients are one. For some neutral aqueous species, such as `CO2(aq)`, predefined models exists and for those
~~~{.cpp}
using namespace Reaktoro; {delete}
ChemicalEditor editor; {delete}
aqueousphase.setActivityModelDrummondCO2();
aqueousphase.setActivityModelSetschenow("NaCl(aq)", 0.1);
~~~

The method `ChemicalEditor::addAqueousPhase` returns an object of class `AqueousPhase`. After its creation, the user can use the `aqueousphase` object to set, for example, the Debye--Huckel model as the chemical model of the aqueous phase, which is used, among other things, to to calculate the activities of the aqueous species:

## Using class Gems {#using-gems}


~~~{.cpp}
using namespace Reaktoro; {delete}
Gems gems("some-gems-project-file.lst");

ChemicalSystem system = editor;
~~~

## Using class Phreeqc {#using-phreeqc}


~~~{.cpp}
using namespace Reaktoro; {delete}
Phreeqc phreeqc;("some-phreeqc-database-file.dat");
phreeqc.execute("some-phreeqc-script-file.dat");

ChemicalSystem system = phreeqc;
~~~

## Using class PhreeqcEditor {#using-phreeqc-editor}

~~~{.cpp}
using namespace Reaktoro; {delete}
PhreeqcEditor editor;("some-phreeqc-database-file.dat");
editor.setAqueousPhase("H O C Na Cl Ca Mg");
editor.setGaseousPhase({"CO2(g)"});
editor.setMineralPhases({"Calcite", "Dolomite"});

ChemicalSystem system = editor;
~~~


[ChemicalEditor]: @ref Reaktoro::ChemicalEditor
[ChemicalSystem]: @ref Reaktoro::ChemicalSystem
