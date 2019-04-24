.. |H2O| replace:: H\ :sub:`2`\ O
.. |H2Og| replace:: H\ :sub:`2`\ O\ (g)
.. |CO2| replace:: CO\ :sub:`2`
.. |CO2g| replace:: CO\ :sub:`2`\ (g)
.. |MgCl2| replace:: MgCl\ :sub:`2`
.. |CaCl2| replace:: CaCl\ :sub:`2`
.. |%vol| replace:: %\ :sub:`vol`
.. |SiO2| replace:: Si0\ :sub:`2`
.. |CaCO3| replace:: CaCO\ :sub:`3`
.. |NaCl| replace:: CaCl
.. |HCl| replace:: HCl
.. |CaMg(CO3)2| replace:: CaMg(CO\ :sub:`3`)\ :sub:`2`

.. |H+| replace:: H\ :sup:`+`
.. |Ca++| replace:: Ca\ :sup:`2+`
.. |Mg++| replace:: Mg\ :sup:`2+`
.. |HCO3-| replace:: HC0\ :sub:`3`\ :sup:`-`
.. |CO2(aq)| replace:: C0\ :sub:`2` (aq)

Kinetic path of calcite reacting with |HCl|-brine
=================================================

In this tutorial, we show how Reaktoro can be used to compute chemical kinetic
calculations with both equilibrium- and kinetically-controlled reactions.
An example below demonstrate this for a simple mineral dissolution modelling,
in which calcite reacts with a 1 molal |HCl|-brine at temperature 298.15 K and
pressure 1 MPa

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step

Below, you find next a step-by-step explanation of the above script.

Importing python packages
-------------------------

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 1
    :end-before: Step 2

Using Reaktoro in Python requires first an import of the python package **reaktoro**.
From this point on, we are able to use the library components of Reaktoro (classes,
methods, constants), which are needed to define our chemical system and chemical
reaction modelling problems.

.. note::

    To simplify the tutorials, we use ``from reaktoro import *``, which imports all
    components of the **reaktoro** package into the default Python namespace. We note
    that this can potentially create name conflicts when used in bigger projects. For
    your applications, consider using ``import reaktoro as rkt`` instead, and call
    classes and methods as ``rkt.Database``, ``rkt.ChemicalSystem``, ``rkt.equilibrate``,
    and etc.

Add phases to the chemical system
---------------------------------

The class ``ChemicalEditor`` is used later to conveniently create instances of classes
``ChemicalSystem`` and ``ReactionSystem``.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 2
    :end-before: Step 3

In particular, we specify **aqueous** and **mineral** phases that should be
considered in the chemical system. The aqueous phase is defined by the mixing
of |H2O|, |HCl|, and |CaCO3|. The mineral phase is defined with only calcite.

.. |supcrt98| replace:: :download:`supcrt98.xml <../../../../databases/supcrt/supcrt98.xml>`
.. |slop98| replace:: :download:`slop98.dat <../../../../databases/supcrt/slop98.dat>`

.. note::
    The default database used by the ``editor`` is the |supcrt98| database.
    It was generated from the original SUPCRT92 database file |slop98|. You are
    welcome to inspect these files and learn more about the chemical species
    available in them. You can also read more about the available thermodynamic
    databases supported in Reaktoro at :ref:`Thermodynamic Databases`.

    All SUPCRT92 thermodynamic databases have been embedded into Reaktoro.
    Thus, you don't actually need to have a database file named ``supcrt98.xml``
    in a local directory when initializing the ``Database`` object. If you want
    to use a customized database file, however, also named ``supcrt98.xml``, then
    your local file will be used instead.

    If you are using a customized version of a thermodynamic database, consider
    changing its name (e.g., ``custom-supcrt98.xml``) to avoid accidental
    use of an embedded database. This can happen if you do not give a correct
    path to your custom database file.

Add reactions to the chemical reactions of the system
-----------------------------------------------------

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 3
    :end-before: Step 4

Next, we define a mineral reaction of the chemical system. In particular, we set
the equation by ``setEquation()`` method. We add the ``MineralMechanism`` to the
``editor`` to prescribe neutral and acidic mechanisms of the mineral reactions.
Finally, we provide the surface area, which appropriate units are checked and
correspondingly set by function ``setSpecificSurfaceArea()``. Here, **logk** is
the kinetic rate constant of the reaction in log scale, whereas **Ea** is the
Arrhenius activation energy. The units of both constants must be provided.

Create chemical and reaction systems
------------------------------------

Create instances of ``ChemicalSystem`` and ``ReactionSystem``.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 4
    :end-before: Step 5

``ChemicalSystem`` is a class that represents a system and its attributes and properties,
such as phases (in our case aqueous and mineral ones), species, elements (5 in total, i.e.,
*H*, *0*, *Ca*, *C*, *Cl*), formula matrix, as well as chemical and thermodynamical model.
Class ``ReactionSystem`` represents a system of the chemical reaction by a collection of
``Reaction`` class instances. It provides convenient methods that calculate the equilibrium
constants, reaction quotients, and rates of the reactions.


Define the partition of the system
----------------------------------

For the partition of a chemical system, we use corresponding class ``Partition``.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 5
    :end-before: Step 6

In this example, the mineral reaction is specified to be under kinetic control by the
``setKineticPhases`` method. The aqueous species are assumed to be in the chemical
equilibrium at all times (capable of reacting instantaneously to a new state of
equilibrium).

.. note::
    In general, the system can be partitioned into *equilibrium*, *kinetic* and
    *inert species*. For the *equilibrium* species, the composition is governed
    by the chemical equilibrium (calculated by the minimization of their Gibbs
    energy subject to some equilibrium constraints). The *kinetic* species are
    the species whose composition is governed by chemical kinetics (by solving a
    system of ordinary differential equations that model the kinetics of a system
    of reactions). The *inert* species are the species whose composition is
    invariable.

In general, the aqueous species react among themselves at much faster rates than
mineral dissolution reactions, and thus this *partial equilibrium assumption* is
plausible and fairly accurate in most cases.


Defining the chemical equilibrium problem
-----------------------------------------

After constructing the chemical system and specifying the partitioning of the
species, we proceed to the definition of the ``EquilibriumProblem``.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 6
    :end-before: Step 7

We specify the above-defined partition of the chemical system by ``setPartition``
method applied to ``EquilibriumProblem``. We prescribe the amount of injected
aqueous fluid resulting from the mixture of 1 kg of water with 1 mole of |HCl|.
The default temperature and the pressure are set to 298.15 K and 1 MPa.

Calculate the chemical equilibrium state
----------------------------------------

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 7
    :end-before: Step 8

In this step, we use the function ``equilibrate`` to calculate the chemical equilibrium
state of the system. For this calculation, Reaktoro uses an efficient **Gibbs energy
minimization** computation to determine the species' amounts that correspond to a
state of minimum Gibbs energy in the system, while satisfying the prescribed amount
conditions for the temperature, pressure, and element amounts. The result is stored in
the objects ``state0`` of class ``ChemicalState``, a computational representation of
the state of a multiphase chemical system defined by its temperature (*T*), pressure
(*P*), and molar composition (*n*).

.. attention::

    In the vast majority of cases, there is one object of class ``ChemicalSystem`` in
    the code and one or more objects of class ``ChemicalState`` describing different
    states of the chemical system. Reaktoro differentiates these two independent
    concepts: *chemical system definition* and *chemical system state*.

Setting the mass of mineral
----------------------------

Prescription of the mass of calcite is performed as follows:

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 8
    :end-before: Step 9

Kinetic path calculations
-------------------------

To be able to run the simulation of the chemical kinetic path, we introduce the
instance of the class ``KineticPath`` that enables this functionality.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 9
    :end-before: Step 10

The instance of kinetic path solver is provided with the partition to the equilibrium,
kinetic, and inert species defined above.

Plotting of evolution of different properties of the chemical system
--------------------------------------------------------------------

Usage of the ``ChemicalPlot`` allows us to create plots illustrating the evolution of
different properties of the chemical system over the time interval.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 10
    :end-before: Step 11

The code above plots four different graphics with evolution of
concentration of calcium, mass of calcite, pH, and concentration of |Ca++| and |HCO3-|.
In particular, instance `plot4` illustrates the depiction of two concentrations with
respect to the time (in hours). Methods ``plot4.x()`` and ``plot4.y()`` fetch properties
that are meant to be plotted on x and y axises, respectively.
The arguments of the function ``plot1.y(""speciesMolality(Ca++ units=mmolal)", "Ca++")``
stand for the name of the property we want to fetch from the chemical state and the tag,
we want to illustrate it with.

Solve the chemical kinetics problem
-----------------------------------

Finally, we solve the kinetic path problem using the method ``solve`` provided the
initial state (here, ``state0``), the initial and final time of the kinetic path, and
time units of the specified time (e.g., `s`, `minute`, `day`, `year`, etc.).

.. literalinclude:: ../../../../demos/python/demo-kineticpath-calcite-hcl.py
    :start-at: Step 11
    :end-before: Step 12


.. _ChemicalState: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html