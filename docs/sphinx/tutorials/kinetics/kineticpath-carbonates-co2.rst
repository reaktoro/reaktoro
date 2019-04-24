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
.. |CaMg(CO3)2| replace:: CaMg(CO\ :sub:`3`)\ :sub:`2`
.. |MgCO3| replace:: MgCO\ :sub:`3`

.. |H+| replace:: H\ :sup:`+`
.. |Ca++| replace:: Ca\ :sup:`2+`
.. |Mg++| replace:: Mg\ :sup:`2+`
.. |HCO3-| replace:: HC0\ :sub:`3`\ :sup:`-`
.. |CO2(aq)| replace:: C0\ :sub:`2` (aq)

Kinetic modelling of |CO2| injection into carbon saline aquifers
================================================================

In this tutorial, we show how Reaktoro can be used to model
water–gas–rock interactions produced by the carbon dioxide injection
by applying both chemical kinetics and chemical equilibrium methodologies.

The injection of carbon dioxide in saline aquifers perturbs
the reservoir and initiates several physical and chemical phenomena
due to the interactions of the injected gas with the resident fluid
and the reservoir rock.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step

You find next a step-by-step explanation of the above script.

Importing python packages
-------------------------

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 1
    :end-before: Step 2

To use Reaktoro in Python, it requires to import the python package
**reaktoro**. From this point on, we are able to use the library components of
Reaktoro (classes, methods, constants), which are needed to define our chemical
system and chemical reaction modelling problems.

.. note::

    To simplify the tutorials, we use ``from reaktoro import *``, which imports
    all components of the **reaktoro** package into the default Python
    namespace. We note that this can potentially create name conflicts when used
    in bigger projects. For your applications, consider using ``import reaktoro as rkt``
    instead, and call classes and methods as ``rkt.Database``, ``rkt.ChemicalSystem``,
    ``rkt.equilibrate``, and etc.

Creating the chemical editor
----------------------------

We start from defining the chemical system using the ``ChemicalEditor`` class:

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 2
    :end-before: Step 3

In particular, we specify **aqueous**, **gaseous**, and **mineral** phases
that should be considered in the chemical system.
The aqueous phase is defined by mixing of |H2O|, |NaCl|, |CaCO3|, and |MgCO3|.
The gaseous phase is represented by |H2Og| and |CO2g|.
The mineral phase contains four mineral species:
calcite |CaCO3|, magnesite |MgCO3|, dolomite |CaMg(CO3)2|, and halite |Nacl|.

.. |supcrt98| replace:: :download:`supcrt98.xml <../../../../databases/supcrt/supcrt98.xml>`
.. |slop98| replace:: :download:`slop98.dat <../../../../databases/supcrt/slop98.dat>`

.. note::
    The default database used by the ``editor`` is the |supcrt98| database.
    It was generated from the original SUPCRT92 database file |slop98|.
    You are welcome to inspect these files and learn more about the
    chemical species available in them. You can also read more about the available
    thermodynamic databases supported in Reaktoro at :ref:`Thermodynamic Databases`.

    All SUPCRT92 thermodynamic databases have been embedded into Reaktoro.
    Thus, you don't actually need to have a database file named ``supcrt98.xml``
    in a local directory when initializing the ``Database``
    object. If you want to use a customized database file, however, also named
    ``supcrt98.xml``, then your local file will be used instead.

Add reactions to the chemical reactions of the system
-----------------------------------------------------

Next, we define minerals' reactions that will be modelled by kinetics
involving calcite, magnesite, and dolomite mineral phases.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 3
    :end-before: Step 4

In particular, we set the equations of the reaction by ``setEquation()`` method. We add
the ``MineralMechanism`` to the ``editor`` to prescribe neutral and acidic mechanisms
of the mineral reactions. Finally, we provide the surface area, which appropriate
units are checked and correspondingly set by function ``setSpecificSurfaceArea()``.
Here, **logk** is the kinetic rate constant of the reaction in log scale, whereas
**Ea** is the Arrhenius activation energy. The units of both constants must be provided.

Create the chemical and reaction systems
----------------------------------------

Create instances of ``ChemicalSystem`` and ``ReactionSystem``.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 4
    :end-before: Step 5

Here, ``ChemicalSystem`` is a class that represents a system as well as its attributes and properties, such as
phases (in our case aqueous and mineral ones), species, elements (in total five elements,
i.e., *H*, *0*, *Ca*, *C*, *Cl*), formula matrix, as well as chemical and thermodynamical model.
Class ``ReactionSystem`` represents a system of the chemical reactions by a collection of ``Reaction`` instances.
It also provides convenient methods that calculate the equilibrium constants, reaction quotients,
and rates of the reactions.


Define the partition of the chemical system
-------------------------------------------

For the partition of a chemical system, we use corresponding class ``Partition``.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 5
    :end-before: Step 6

In this example, three mineral reactions are specified to be under kinetic control by the
``setKineticPhases`` method. The aqueous and gaseous species, as well as halite, are assumed
to be in the chemical equilibrium at all times (capable of reacting instantaneously
to a new state of equilibrium).

.. note::
    In general, the system can be partitioned into *equilibrium*, *kinetic* and *inert species*.
    For the *equilibrium* species, the composition is governed by the chemical equilibrium
    (calculated by the minimization of their Gibbs energy subject to some equilibrium constraints).
    The *kinetic* species are the species whose composition is governed by chemical kinetics
    (by solving a system of ordinary differential equations that models the kinetics of a system of reactions).
    The *inert* species are the species whose composition is invariable.

Since the aqueous species react among themselves at much faster rates than
mineral dissolution reactions, and thus this *partial equilibrium assumption* is
plausible, and fairly accurate in most cases.

Defining the chemical equilibrium problem
-----------------------------------------

After constructing the chemical system and specifying the partitioning of the species,
we can proceed to the *equilibrium problem*.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 6
    :end-before: Step 7

The above-defined partition of the chemical is set by ``setPartition()`` method.
We prescribe the amount of injected aqueous fluid resulting from
the mixture of 1 kg of water with 1 mole of |NaCl|. The amount of injected |CO2| is set to be
1 mol. The default temperature and the pressure are set to 298.15 K and 1 MPa.

Calculate the chemical equilibrium state
-----------------------------------------

To calculate the equilibrium state of the system comprised of the injected carbon
dioxide, subsurface fluid, and the rock-forming minerals, one uses the following method:

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 7
    :end-before: Step 8

For this calculation, Reaktoro uses an efficient
**Gibbs energy minimization** computation to determine the species' amounts that
correspond to a state of minimum Gibbs energy in the system, while satisfying
the prescribed amount conditions for the temperature, pressure, and element
amounts. The result is stored in the objects ``state0`` of class ``ChemicalState``,
a computational representation of the state of a multiphase chemical system defined
by its temperature (*T*), pressure (*P*), and molar composition (*n*).

.. attention::

    In the vast majority of cases, there is one object of class ``ChemicalSystem``
    in the code and one or more objects of class ``ChemicalState`` describing different
    states of the chemical system. Reaktoro differentiates these two independent concepts:
    *chemical system definition* and *chemical system state*.

Setting the mass of mineral
----------------------------

Prescription of the mass of the rock-forming minerals that are in disequilibrium with the
|CO2| saturated subsurface fluid.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 8
    :end-before: Step 9

Setting the  kinetic path calculations
--------------------------------------

To be able to run the simulation of a kinetic process of mineral dissolution/precipitation,
we introduce the instance of the class ``KineticPath`` that enables this functionality.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 9
    :end-before: Step 10

The instance of kinetic path solver is provided with the partition to the equilibrium,
kinetic, and inert species defined above.

Plotting of evolution of different properties
---------------------------------------------

Usage of the ``ChemicalPlot`` allows us to create plots illustrating the evolution of
different properties of the chemical system over the time interval.

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 10
    :end-before: Step 11

The code above plots four different graphics with the evolution of pH, the concentration
of calcium and magnesium, and the mass of calcite dolomite. In particular, instance
`plot1` illustrates the depiction of two properties, i.e., calcium and magnesium
concentrations, with respect to the time (in hours). Methods ``plot1.x()`` and
``plot1.y()`` fetch properties that are meant to be plotted on x and y axises, respectively.
The arguments of the function ``plot1.y("elementMolality(Ca)", "Ca")`` stand for the
name of the property we want to fetch from the chemical state and the tag, we want to
illustrate it with.


Solve the chemical kinetics problem
-----------------------------------

Finally, we calculate the transient state of the entire system compromised of the
rock-forming minerals, the subsurface fluid, and the emerged |CO2|-rich phase. This
is performed by using the method ``solve`` provided the
initial state (here, ``state0``), the initial and final time of the kinetic path, and
time units of the specified time (e.g., `s`, `minute`, `hour`, `day`, `year`, etc.).

.. literalinclude:: ../../../../demos/python/demo-kineticpath-carbonates-co2.py
    :start-at: Step 11
    :end-before: Step 12

.. _ChemicalState: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html