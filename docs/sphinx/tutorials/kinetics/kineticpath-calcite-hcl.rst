.. |degC| replace:: °C
.. |H2O| replace:: H\ :sub:`2`\ O
.. |H2Og| replace:: H\ :sub:`2`\ O\ (g)
.. |CO2| replace:: CO\ :sub:`2`
.. |CO2g| replace:: CO\ :sub:`2`\ (g)
.. |MgCl2| replace:: MgCl\ :sub:`2`
.. |CaCl2| replace:: CaCl\ :sub:`2`
.. |%vol| replace:: %\ :sub:`vol`
.. |SiO2| replace:: Si0\ :sub:`2`
.. |CaCO3| replace:: CaCO\ :sub:`3`
.. |NaCl| replace:: NaCl
.. |HCl| replace:: HCl
.. |CaMg(CO3)2| replace:: CaMg(CO\ :sub:`3`)\ :sub:`2`
.. |H+| replace:: H\ :sup:`+`
.. |Ca++| replace:: Ca\ :sup:`2+`
.. |Mg++| replace:: Mg\ :sup:`2+`
.. |HCO3-| replace:: HC0\ :sub:`3`\ :sup:`-`
.. |CO2(aq)| replace:: C0\ :sub:`2` (aq)

Dissolution of calcite in an acidic |HCl|-solution
==================================================

This tutorial demonstrates how Reaktoro can be used for modelling the dissolution
of calcite in an acidic |HCl|-solution at temperature 30 |degC| and pressure 1
bar using chemical kinetics. A partial equilibrium assumption is considered
here so that aqueous species react using a chemical equilibrium model, while
calcite reacts with the aqueous solution using a chemical kinetics model.

The full script for the calculation is shown below, followed by a step-by-step
explanation afterwards.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step

Importing reaktoro Python package
---------------------------------

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 1
    :end-before: Step 2

Here, we again import the **reaktoro** Python package so that we can use its
classes and methods for performing the chemical reaction calculations.

Defining the phases in the chemical system
------------------------------------------

The class `ChemicalEditor`_ is used to conveniently create instances of
classes `ChemicalSystem`_ and `ReactionSystem`_.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 2
    :end-before: Step 3

In particular, we specify **aqueous** and **mineral** phases that should be
considered in the chemical system. The aqueous phase is defined by the mixing
of |H2O|, |HCl|, and |CaCO3| (effectively, collecting all aqueous species in the
database that contains elements H, O, C, Cl, and Ca, which are the elements in this
list of compounds). There is only one pure mineral phase: the calcite phase.

Defining kinetically-controlled reactions
-----------------------------------------

Next, we define a mineral reaction for calcite and set its rate parameters.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 3
    :end-before: Step 4

We set the reaction equation using the ``setEquation`` method of the
class `MineralReaction`_.

.. todo::
    The current need for setting a reaction equation for mineral reactions will
    be removed in the future, which will simplify the setup of the problem.

Then we add two mineral kinetic mechanisms for the reaction: neutral and
acidic. This is done with the ``addMechanism`` method, where we set, for
example, ``logk``, the kinetic rate constant of the reaction in log scale, and
``Ea``, the Arrhenius activation energy. The values shown for ``logk`` and
``Ea`` were collected from:

*Palandri, J.L., Kharaka, Y.K. (2004). A compilation of rate parameters of
water-mineral interaction kinetics for application to geochemical modeling.
U.S. Geological Survey Open File Report (Vol. 2004–1068). Menlo Park,
California*.

.. note::
    The units of ``logk`` and ``Ea`` can be different, but
    convertible to ``mol/(m2*s)`` and ``J/mol`` respectively.

Finally, we provide the specific surface area of the mineral using method
``setSpecificSurfaceArea`` of class `MineralReaction`_, which can be specified
in units of ``m2/g`` or ``m2/m3``. Compatible units are allowed, such as
``cm2/mg`` or ``m2/dm3``, and combinations.

Creating the chemical and reaction systems
------------------------------------------

Create instances of `ChemicalSystem`_ and `ReactionSystem`_.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 4
    :end-before: Step 5

`ChemicalSystem`_ is a class that represents a system and its attributes and properties,
such as phases (in our case aqueous and mineral ones), species, elements (5 in total, i.e.,
*H*, *O*, *Ca*, *C*, *Cl*), formula matrix, as well as chemical and thermodynamical model.
Class `ReactionSystem`_ serves as a system of the chemical reaction by a collection of
`Reaction`_ class instances. It provides convenient methods that calculate the equilibrium
constants, reaction quotients, and rates of the reactions.


Specifying the equilibrium and kinetic species
----------------------------------------------

For the partition of a chemical system into equilibrium and kinetic species, we
use the class `Partition`_. We only need to specify which species are kinetic
species, and all others will be equilibrium species by default.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 5
    :end-before: Step 6

We set species ``Calcite`` (the only species in the mineral phase also called
``Calcite``!) to be the only kinetic species. This will allow us to model the
dissolution of calcite using chemical kinetics, while all other species (the
aqueous species) are modelled using chemical equilibrium (i.e., their amounts
are updated over time using chemical equilibrium calculations).

.. note::
    In general, the chemical system can be partitioned into *equilibrium*,
    *kinetic* and *inert species*. For the *equilibrium species*, their amounts
    are governed by chemical equilibrium (calculated by solving a Gibbs energy
    minimization problem). The *kinetic species* are the species whose amounts
    are governed by chemical kinetics (calculated by solving a system of
    ordinary differential equations that model the kinetics of a system of
    reactions). The *inert species* are the species whose amounts are constant.

.. attention::
    Aqueous species react among themselves at much faster rates than mineral
    dissolution reactions in general, and thus this *partial equilibrium
    assumption* is plausible and reasonably accurate.

Defining the chemical equilibrium problem
-----------------------------------------

After constructing the chemical system and specifying the partitioning of the
species, we proceed to the definition of the `EquilibriumProblem`_.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 6
    :end-before: Step 7

We specify the equilibrium/kinetic partitioning of the chemical system using
method ``setPartition`` of class `EquilibriumProblem`_. We then prescribe what
should be the initial state of the equilibrium species (the aqueous species in
this case), before we start the chemical kinetics calculation that will
simulate the dissolution of calcite in this aqueous fluid.

By mixing 1 kg of |H2O| and 1 mmol of |HCl| at 30 |degC| and 1 bar, we should
produce a chemical equilibrium state that corresponds to an acidic aqueous
fluid. The species in this fluid will be in disequilibrium with ``Calcite``
(our single kinetic species in this setup) since only equilibrium species
(i.e., the aqueous species) are considered during the next chemical
equilibrium calculation.

Calculating the initial chemical equilibrium state of the fluid
---------------------------------------------------------------

We now use the function ``equilibrate`` to calculate the chemical equilibrium
state of the equilibrium partition, not the entire chemical system.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 7
    :end-before: Step 8

For this calculation, Reaktoro uses an efficient **Gibbs energy minimization**
algorithm to determine the amounts of the equilibrium species that correspond
to a state of minimum Gibbs energy in the equilibrium partition only, at given
conditions of temperature, pressure, and element amounts in the equilibrium
partition. The result is stored in the object ``state0`` of class
`ChemicalState`_, a computational representation of the state of a multiphase
chemical system defined by its temperature (*T*), pressure (*P*), and vector of
species amounts (*n*).

.. note::
    Consider using class `EquilibriumSolver`_ instead of method ``equilibrate``
    when a sequence of chemical equilibrium calculations are needed for
    increased performance. The method ``equilibrate`` is a convenient function
    for solving chemical equilibrium problems, but has some overhead (e.g., an
    object of class `EquilibriumSolver`_ is created in each call).

Setting the initial mass of calcite
-----------------------------------

To simulate the kinetic dissolution of calcite in the aqueous fluid we defined
before, we need to specify its initial amount. Below, we set the initial mass
of species ``Calcite`` to 100 g.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 8
    :end-before: Step 9

Performing the kinetic path calculation
---------------------------------------

To be able to run the simulation of the chemical kinetic path, we use class
`KineticPath`_.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 9
    :end-before: Step 10

Note that here again, we need to specify the partitioning of the chemical system
into equilibrium, kinetic, and inert species.

.. attention::
    For repeated chemical kinetics calculations (e.g., in a reactive transport
    simulation, where kinetics is performed for each mesh cell/node), consider
    using the class `KineticSolver`_ instead for avoiding some overhead of
    class `KineticPath`_.

Plotting the evolution of different properties of the chemical system
---------------------------------------------------------------------

Usage of the `ChemicalPlot`_ allows us to create plots that show the evolution of
different properties of the chemical system during the reactive process.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 10
    :end-before: Step 11

The code above plots four different graphics that show the concentration of
calcium, mass of calcite, pH, and concentration of ions |Ca++| and |HCO3-| over
time. In particular, ``plot4`` shows the concentrations of ions |Ca++| and
|HCO3-| with respect to the time (in hours). Methods ``plot4.x()`` and
``plot4.y()`` fetch properties that are meant to be plotted on x and y axes,
respectively. The arguments of the function ``plot1.y(""speciesMolality(Ca++
units=mmolal)", "Ca++")`` stand for the name of the property from the chemical
state we want to plot and the label to be used in the plot.

.. note::
    A list of all possible quantities that can be plotted is shown in the class
    `ChemicalQuantity`_, which provides an interface for convenient ways of their
    retrieval.

Solving the chemical kinetics problem
-------------------------------------

Finally, we solve the kinetic path problem.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 11
    :end-before: Step 12

This step executes the method ``solve`` of class `KineticPath`_, which
requires the initial state of the system (100 g of calcite in disequilibrium
with a 1 mmolal HCl aqueous solution at 30 |degC| and 1 bar, represented with
the object ``state0``), the initial and final time of the kinetic path
calculation (``t0`` and ``t1``, respectively), and the time unit of the
specified time parameters (e.g., `s`, `minute`, `day`, `year`, etc.).


.. _ChemicalEditor: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html
.. _ChemicalPlot: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalPlot.html
.. _ChemicalQuantity: https://www.reaktoro.org/cpp/classReaktoro_1_1ChemicalQuantity.html
.. _ChemicalState: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html
.. _ChemicalSystem: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html
.. _EquilibriumProblem: https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumProblem.html
.. _EquilibriumSolver: https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html
.. _KineticPath: https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html
.. _KineticSolver: https://reaktoro.org/cpp/classReaktoro_1_1KineticSolver.html
.. _MineralMechanism: https://reaktoro.org/cpp/structReaktoro_1_1MineralMechanism.html
.. _MineralReaction: https://reaktoro.org/cpp/classReaktoro_1_1MineralReaction.html
.. _Partition: https://reaktoro.org/cpp/classReaktoro_1_1Partition.html
.. _Reaction: https://reaktoro.org/cpp/classReaktoro_1_1Reaction.html
.. _ReactionSystem: https://reaktoro.org/cpp/classReaktoro_1_1ReactionSystem.html