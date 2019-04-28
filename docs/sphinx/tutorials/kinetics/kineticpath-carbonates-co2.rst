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
.. |MgCO3| replace:: MgCO\ :sub:`3`
.. |CaMg(CO3)2| replace:: CaMg(CO\ :sub:`3`)\ :sub:`2`
.. |H+| replace:: H\ :sup:`+`
.. |Ca++| replace:: Ca\ :sup:`2+`
.. |Mg++| replace:: Mg\ :sup:`2+`
.. |HCO3-| replace:: HC0\ :sub:`3`\ :sup:`-`
.. |CO2(aq)| replace:: C0\ :sub:`2` (aq)

Dissolution of carbonate minerals in a |CO2|-saturated brine
============================================================

In this tutorial, we demonstrate how Reaktoro can be used to kinetically model
the dissolution of carbonate minerals (calcite, magnesite and dolomite) in a
|CO2|-saturated brine.

Below is the entire Python script showing the calculation, followed by a
step-by-step explanation. You may want to just have a quick look at this
relatively long script, and then proceed to the next sections for a more
detailed discussion!

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step

Importing the reaktoro Python package
-------------------------------------

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 1
    :end-before: Step 2

Here we import the **reaktoro** Python package to enable us to use all its library
components (classes, methods, constants).

Defining the chemical system
----------------------------

In this simulation, we consider an *aqueous phase* (to model the brine
solution), a *gaseous phase* (to model the |CO2|-rich phase with water vapour),
and four pure *mineral phases*: halite, calcite, magnesite, and dolomite. These
are the phases that will either exist initially during the simulation, or that
could potentially appear later as the calculations proceed.

.. attention::
    All potential phases that could appear in a reactive process should ideally
    be considered when defining the chemical system. If one or more of these
    phases are ignored, then the equilibrium and kinetics calculations cannot
    identify if they should actually exist or not with positive amounts.
    Unrealistic results may occur, such as, for example, an aqueous phase
    containing more |CO2| dissolved than it could, because a gaseous phase,
    which should contain the excess of |CO2|, was not considered in the
    chemical system.

The code below uses class `ChemicalEditor`_ to define our chemical system with
the phases of interest and their species:

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 2
    :end-before: Step 3

The aqueous phase is defined by considering all aqueous species in the database
that could form once the substances |H2O|, |NaCl|, |CaCO3|, and |MgCO3| are
mixed. The gaseous phase is defined so that only the gaseous species |H2Og| and
|CO2g| are considered. There are four pure mineral phases: calcite (|CaCO3|),
magnesite (|MgCO3|), dolomite (|CaMg(CO3)2|), and halite (|Nacl|).

Defining the kinetically-controlled reactions
---------------------------------------------

A *partial equilibrium assumption* is considered in this simulation. This
simplifies the problem by assuming that those species that react at much faster
rates are *continually in chemical equilibrium*. They are referred to as
*equilibrium species*. The remaining species (those reacting at relatively
slower rates) are referred to as *kinetic species*.

Because aqueous and gaseous species, as well as halite, react at relatively
fast rates, they are reasonable candidates in this problem for being equilibrium
species. The kinetic species are thus the carbonate minerals: calcite, magnesite,
and dolomite.

With this partial equilibrium assumption, there is no need to specify kinetic
rate laws for the fast reacting equilibrium species. For these species,
chemical equilibrium calculations are performed to update their amounts as the
amounts of each kinetic species change over time (i.e., the equilibrium species
react instantaneously to a new state of equilibrium as the kinetic species
react and perturb the current partial equilibrium state).

Thus, we only need to define kinetic rates for the relatively slow reacting
carbonate minerals (our kinetic species in this simulation), which is shown
below:

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 3
    :end-before: Step 4

We set the equation of each mineral reaction using method ``setEquation`` of
class `MineralReaction`_.

.. todo::
    The current need for setting a reaction equation for mineral reactions will
    be removed in the future, which will simplify the setup of the problem.

We then prescribe neutral and acidic mechanisms for each mineral reaction using
the method ``addMechanism`` of class `MineralMechanism`_, using values for
``logk``, the kinetic rate constant of the reaction in log scale, and ``Ea``,
the Arrhenius activation energy. The units of both parameters must be provided
as shown in the example, and other compatible units are allowed.

Finally, we define the specific surface area of the mineral using method
``setSpecificSurfaceArea`` of class `MineralReaction`_. Any units compatible to
``m2/kg`` or ``m2/m3`` are allowed (e.g., ``cm2/g``, ``mm2/mm3``).

.. note::
    The values shown for ``logk`` and ``Ea`` were collected from

    *Palandri,
    J.L., Kharaka, Y.K. (2004). A compilation of rate parameters of
    water-mineral interaction kinetics for application to geochemical modeling.
    U.S. Geological Survey Open File Report (Vol. 2004–1068). Menlo Park,
    California*.

Creating the chemical and reaction systems
------------------------------------------

Create instances of `ChemicalSystem`_ and `ReactionSystem`_.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 4
    :end-before: Step 5

Class `ChemicalSystem`_ is used to represent a chemical system, containing one
or more phases  (in this problem, aqueous, gaseous, and mineral phases), each
phase containing one or more species (the aqueous phase containing many aqueous
species, the gaseous phase containing two, and each mineral phase containing a
single mineral species with the same name as the name of the phase). This class
is also used to calculate thermodynamic properties of the phases and species
(standard thermodynamic properties, activities, chemical potentials, phase
molar volume, etc.).

Class `ReactionSystem`_ is a collection of `Reaction`_ objects used to
represent a system of chemical reactions that are controlled by chemical
kinetics. These classes provide convenient methods for the calculation of
equilibrium constants, reaction quotients, and rates of the reactions.

Specifying equilibrium and kinetic species
------------------------------------------

In Reaktoro, the species in a chemical system can be partitioned into groups of:

* *equilibrium species*,
* *kinetic species*, and
* *inert species*.

For the *equilibrium species*, their amounts are governed by chemical
equilibrium (calculated via a Gibbs energy minimization). The amount of *kinetic
species* are controlled by chemical kinetics (by solving a system
of ordinary differential equations that models the kinetics of a system of
reactions). The *inert species* maintain their amounts constant in all chemical
equilibrium or kinetics calculations.

This classification of species can be done using class `Partition`_. By
default, all species are considered to be equilibrium species in this class.
Thus, we only need to specify which ones are kinetic species:

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 5
    :end-before: Step 6

In this case, the mineral species calcite, magnesite and dolomite are specified
to be *kinetic species* using method ``setKineticSpecies`` of class `Partition`_.

.. note::
    Method ``setKineticPhases`` of class `Partition`_ could also be used here.
    This method sets all species in the given phases to be kinetic species, and
    it is more convenient if a phase has many species. However, since each of
    the mineral phases considered here only contains a single mineral species,
    method ``setKineticSpecies`` is a convenient alternative.

Defining the initial state of the equilibrium species
-----------------------------------------------------

In a chemical kinetics calculation, an *initial condition* is needed for the
amounts of both equilibrium and kinetic species.

The equilibrium problem formulated below, using class `EquilibriumProblem`_, is
done so that the initial condition for the amounts of each equilibrium species
result from the solution of a chemical equilibrium problem in which 1 kg of
water is mixed with 1 mol of |NaCl| and 1 mol of |CO2| at 60 |degC| and 100
bar. This amount of |CO2| is sufficient to saturate the brine solution. The
excess will exist in the |CO2|-rich gaseous phase.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 6
    :end-before: Step 7

.. attention::
    To ensure that that the equilibrium calculation performed in the next step
    ignores the kinetic species in the system so that we maintain a
    disequilibrium state between equilibrium and kinetic species, it is
    important not to forget the following command:

    .. code-block:: python

        problem.setPartition(partition)

    Ignoring this step will produce an initial condition for the amounts of
    equilibrium and kinetic species that correspond to a complete equilibrium
    state in the system so that no kinetics calculation makes sense
    afterwards!

Calculating the initial amounts of the equilibrium species
----------------------------------------------------------

We now use the convenient ``equilibrate`` function to calculate the amounts of
the equilibrium species by minimizing the Gibbs energy of the equilibrium
partition only, and not of the entire system. The result is stored in the
object ``state0`` of class `ChemicalState`_, a computational representation of
the state of a multiphase chemical system defined by its temperature (*T*),
pressure (*P*), and vector of species amounts (*n*). We then output this
chemical state to a file.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 7
    :end-before: Step 8

Setting the initial mass of minerals
------------------------------------

We have now to prescribe the initial amounts of the kinetic species (i.e., the
carbonate minerals). This is done below by setting the initial mass of calcite
to 100 g and of dolomite to 50 g. The initial amount of magnesite is zero.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 8
    :end-before: Step 9

Setting the kinetic path calculations
-------------------------------------

To be able to run the simulation of the kinetic process of mineral
dissolution/precipitation, we introduce the instance of the class
`KineticPath`_ that enables this functionality.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 9
    :end-before: Step 10

The instance of kinetic path solver is provided with the partition to the
equilibrium, kinetic, and inert species defined above.

.. attention::
    For repeated chemical kinetics calculations (e.g., in a reactive transport
    simulation, where kinetics is performed for each mesh cell/node), consider
    using the class `KineticSolver`_ instead for avoiding some overhead of
    class `KineticPath`_. Consider also the class `EquilibriumSolver`_,
    instead of method ``equilibrate``, for similar reasons.

Plotting the evolution of thermochemical properties
---------------------------------------------------

Usage of the `ChemicalPlot`_ allows us to create plots illustrating the
evolution of different properties of the chemical system over the time interval
in which kinetics is calculated.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 10
    :end-before: Step 11

The code above plots four different figures with the evolution of pH, the
concentration of calcium and magnesium, and the mass of calcite dolomite. In
particular, `plot1` shows how two properties, i.e., calcium and magnesium
concentrations, with respect to the time (in hours). Methods ``plot1.x()`` and
``plot1.y()`` fetch properties that are meant to be plotted on x and y axises,
respectively. The arguments of the function ``plot1.y("elementMolality(Ca)",
"Ca")`` stand for the name of the property we want to fetch from the chemical
state and the tag, we want to illustrate it with.

.. note::
    A list of all possible quantities that can be plotted is shown in the class
    `ChemicalQuantity`_.

Solving the chemical kinetics problem
-------------------------------------

Finally, we calculate the transient state of the entire system comprised of the
carbonate minerals, the |CO2|-saturated brine fluid, and the |CO2|-rich gaseous
phase. This is performed by using the method ``solve`` of the `KineticPath`_
class, which expects the initial state (``state0``), the initial and final
times (``t0`` and ``t1`` respectively) of the kinetic path, and time units of
the specified time (e.g., `s`, `minute`, `hour`, `day`, `year`, etc.).

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 11
    :end-before: Step 12

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