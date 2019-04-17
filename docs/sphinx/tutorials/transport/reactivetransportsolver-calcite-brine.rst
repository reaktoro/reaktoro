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

.. |H+| replace:: H\ :sup:`+`
.. |Ca++| replace:: Ca\ :sup:`2+`
.. |Mg++| replace:: Mg\ :sup:`2+`
.. |HCO3-| replace:: HC0\ :sub:`3`\ :sup:`-`
.. |CO2(aq)| replace:: C0\ :sub:`2` (aq)

.. |slop98| replace:: :download:`slop98.dat <../../../../databases/supcrt/slop98.dat>`

Reactive transport modelling of a rock core after brine injection (using  ReactiveTrasportSolver class)
=======================================================================================================

In this tutorial, we show how Reaktoro can be used for sequential calculations of the
reactive transport of a rock core after injecting the fluid and rock composition. To do that,
we exploit a predefined class ``ReactiveTransportSolver`` that provides the functionality of
solving advection-diffusion equation on the each step of the time-stepping cycle and logs the
results (customised by the user) into the sequence of txt file (corresponding to each time-step).

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step

You find next a step-by-step explanation of the above script.

Importing the reaktoro package
------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 1
    :end-before: Step 2

Using Reaktoro in Python requires first an import of the python package
**reaktoro**. From this point on, we are able to use the library components of
Reaktoro (classes, methods, constants), which are needed to define our chemical
system and chemical reaction modeling problems.

.. note::

    To simplify the tutorials, we use ``from reaktoro import *``, which imports
    all components of the **reaktoro** package into the default Python
    namespace. We note that this can potentially create name conflicts when used
    in bigger projects. For your applications, consider using ``import reaktoro as rkt``
    instead, and call classes and methods as ``rkt.Database``, ``rkt.ChemicalSystem``,
    ``rkt.equilibrate``, and etc.

Initializing an auxiliary time-related constants
------------------------------------------------

In this step:

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 2
    :end-before: Step 3

we initialize an auxiliary time-related constants from seconds up to years.

Define parameters for the reactive transport simulation
-------------------------------------------------------

Next, we define reactive transport and numerical discretization parameters.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 3
    :end-before: Step 4

First, we define the considered rock domain by setting the coordinate of its left boundary to
0.0 and right boundary to 100.0. The discretization parameters, i.e., both the number of cells and
steps in time, are set to 100. The reactive transport modeling procedure assumes a constant
fluid velocity of v = 1 m/day (1.16 · 10−5 m/s) and the same diffusion coefficient D = 10−9 m2/s
for all fluid species, without dispersivity. The size of the time-step is set to be half a day.
The temperature and pressure of the fluids are 60 °C and 100 bar, respectively.

Defining the chemical system using chemical editor class
--------------------------------------------------------

Reaktoro is a general-purpose chemical solver that avoids as much as possible
presuming specific assumptions about your problems. Thus, you need to specify
how your chemical system should be defined. This encompasses the specification of
all phases in the system as well as the chemical species that compose each
phase. By using the ``ChemicalEditor`` class, you can conveniently achieve
this as shown below:

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 4
    :end-before: Step 5

In this step, we create an object of class ``ChemicalEditor`` and specify two phases,
an **aqueous** and a **mineral** one, should be considered in the chemical system.
The aqueous phase is defined by using a list of compounds, which will be broken
into a list of element names, and the database will then be searched for all species that
could be formed out of those elements. The mineral phase is defined by three mineral species:
quartz |SiO2|, calcite |CaCO3|, and dolomite |CaMg(CO3)2|.

.. note::

    An automatic search for chemical species can result in a large number of
    species in the phase, potentially causing the chemical reaction
    calculations to be more computationally expensive. If you are using
    Reaktoro in demanding applications (e.g., as a chemical solver in a
    reactive transport simulator), you might want to manually specify the
    chemical species of each phase in your chemical system. This can be
    achieved by providing a list of species names as shown below:

    .. code-block:: python

        editor.addAqueousPhaseWithElements([
            'H2O(l)',
            'H+',
            'OH-',
            'Na+',
            'Cl-',
            'Ca++',
            'Mg++',
            'HCO3-',
            'CO2(aq)',
            'CO3--'
            ])

Constructing the chemical system
--------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 5
    :end-before: Step 6

This step is where we create an object of class ``ChemicalSystem`` using the
chemical system definition details stored in the object ``editor``.

.. note::

    ``ChemicalSystem`` is perhaps the main class in Reaktoro. An object of this
    class stores the phases, species and elements in our defined chemical
    system, as well as provides the means to compute many types of thermodynamic
    properties, such as *standard thermodynamic properties* (e.g., standard
    Gibbs energies, standard enthalpies, and standard volumes of species), and
    *thermochemical properties* (e.g., activity and activity coefficients of
    species; density, enthalpy and internal energy of phases). As you learn
    more about other Reaktoro's classes, you will note that an object of class
    ``ChemicalSystem`` is almost always needed for their initialization!

.. note::

    The *activity coefficients* of the aqueous species are calculated using the
    *HKF extended Debye-Hückel model* for solvent water and ionic species, except
    for the aqueous species |CO2| (aq), for which the *Drummond model* is used.

    The *standard chemical potentials* of the species are calculated using the equations
    of state of
    Helgeson and Kirkham,
    Helgeson et al.,
    Tanger and Helgeson [51],
    Shock and Helgeson [52] and
    Shock et al. [53].

    The database file |slop98| from the software SUPCRT92 [54] is used to obtain
    the parameters for the equations of state.

    The equation of state of Wagner and Pruss [55] is used to calculate the
    *density of water* and its temperature and pressure derivatives. Kinetics of
    *dissolution* and *precipitation* of both calcite and dolomite is neglected, i.e.,
    the local equilibrium assumption is employed.


Defining the initial condition of the reactive transport modeling problem
-------------------------------------------------------------------------

We have now defined and constructed our chemical system of interest, enabling
us to move on to the next step in Reaktoro's modeling workflow: *defining our
chemical reaction problems*. Below we define its **initial condition** with already
prescribed equilibrium conditions for *temperature*, *pressure*, and *amounts
of elements* that are consistent with an intention of modeling reactive transport
of injected |NaCl|-|MgCl2|-|CaCl2| brine into the rock-fluid composition of quartz
and calcite at 60 °C and 100 bar. In particular, we consider resident fluid is a 0.7
molal |NaCl| brine in equilibrium with the rock minerals with a calculated pH of 10.0.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 6
    :end-before: Step 7

.. note::
    Our **0.7 molal NaCl brine** is here represented by the mixture of
    1 kg of |H2O| and 0.7 mol of |NaCl|.

.. note::
    After providing the amounts of substances
    |H2O|, |NaCl|, quartz |SiO2|, and calcite |CaCO3| in the above code, Reaktoro parses
    these chemical formulas and determines the elements and their coefficients.
    Once this is done, the amount of each element stored inside the object
    ``problem_ic`` is incremented according to the given amount of substance and
    its coefficient in the formula. The amounts of elements you provide are
    then used as constraints for the Gibbs energy minimization calculation when
    computing the state of chemical equilibrium (i.e., when we try to find the
    amounts of all species in the system that corresponds to a state of minimum
    Gibbs energy and at the same time satisfying the *element amounts
    constraints*).

.. note::
    Please note that we are not condemning the input form shown above in terms
    of element amounts, but only telling you to be attentive with the values
    you input. If you are using Reaktoro as a chemical reaction solver in a
    reactive transport simulator, for example, you'll most likely need to work
    directly with given amounts of elements, which shows that this input form
    is required in certain cases. For such time-dependent modeling problems, you
    often only need to ensure that the initial conditions for elements amounts
    result in feasible initial species amounts.


Defining the boundary condition of the reactive transport modeling problem
--------------------------------------------------------------------------

Next, we define the **boundary condition** of the constructed chemical system with
its *temperature*, *pressure*, and *amounts of elements*.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 7
    :end-before: Step 8

In particular, we prescribe the amount of injected aqueous fluid resulting from
the mixture of
1 kg of water with
0.90 moles of |NaCl|,
0.05 moles of |MgCl2|,
0.01 moles of |CaCl2|, and
0.75 moles of |CO2|, in a state very close to |CO2| saturation.
The temperature and the pressure stay the same, i.e., 60 °C and 100 bar,
respectively.

Calculate the equilibrium states for the initial and boundary conditions
------------------------------------------------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 8
    :end-before: Step 9

In this step, we use the ``equilibrate`` function to calculate the chemical
equilibrium state of the system with the given initial and boundary equilibrium
conditions stored in the object ``problem_ic`` and ``problem_bc``.
For this calculation, Reaktoro uses an efficient
**Gibbs energy minimization** computation to determine the species amounts that
correspond to a state of minimum Gibbs energy in the system, while satisfying
the prescribed amount conditions for temperature, pressure, and element
amounts. The result is stored in the objects ``state_ic`` and ``state_bc`` of
class ``ChemicalState``.

.. attention::

    In the vast majority of cases, there is one object of class ``ChemicalSystem``
    in the code and one or more objects of class ``ChemicalState`` describing different
    states of the chemical system. Reaktoro differentiates these two independent concepts:
    *chemical system definition* and *chemical system state*.


Scaling the phases in the initial condition as required
-------------------------------------------------------

Here, we scale the phases in the initial condition according to the following composition:
98|%vol| |SiO2| (quartz) and 2|%vol| |CaCO3| (calcite) with the porosity of 10%.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 9
    :end-before: Step 10


Scaling the boundary condition state to 1 m3
--------------------------------------------

Next, we scale the boundary condition state to 1 m3

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 10
    :end-before: Step 11

Creating the mesh for the column
--------------------------------

To define the mesh a class ``Mesh`` to initialize it for a class ``ReactiveTransportSolver``.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 11
    :end-before: Step 12

The class accepts the number of cells on the computational domain as well as x-coordinates
of the left and right boundaries (in m). By default, the number of cells is set to 10, whereas
the domain is set to the unit interval.

Creating a chemical field object
--------------------------------

For initializing the reactive transport modelling class, we need to define an object
of a class ``ChemicalField`` with every cell having a state given by state_ic.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 12
    :end-before: Step 13


.. note::
   Alternatively, the chemical field can be initialized by the chemical system common to all
   degrees of freedom in the chemical field.


Defining the reactive transport modeling
----------------------------------------

Finally, we defined the object responsible for the solving reactive transport problem, which
is handled by the class ``ReactiveTransportSolver``.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 13
    :end-before: Step 14

The object is initialized by the chemical system common to all degrees of freedom in the chemical
field. We also provide it with other discretization parameters such as mesh, velocity, diffusion
coefficient, state on the boundary condition, and size of the step for incremental time stepping.
Lastly, we initialize the reactive solver object with the chemical field object specified on the
previous step.

Define the output quantities
----------------------------

Before running time-dependent simulations, we defined another object provided by the class
``ChemicalOutput`` to output the state for every cell, every time step.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 14
    :end-before: Step 15

We provide the name of the file ``reativetransport.txt``, where the data should be collected.
We all have to specify the parameters that we are interested to output. In this case, it is
pH, molality of |H+|, |Ca++|, |Mg++|, |HCO3-|, |CO2(aq)|, as well as phase volume of calcite
and dolomite.


Running the reactive transport simulations in the cycle
-------------------------------------------------------

Before proceeding to the simulation of reactive transport in the considered interval,
we set the initial time and a counter for the step considered in this cycle.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 15
    :end-before: Step 16

The cycle for the reactive transport simulation proceeds until we haven't made all the steps
in time. At each time step, we print the progress of the simulations, which are performed by
the class ``ReactiveTransportSolver``. Each call of function ``rt.step`` performs one reactive
transport time-step, i.e., solving of the advection-diffusion problem using ``ReactiveTransport``
class and writing the results in the file ``reativetransport-step.txt``, where ``step`` indicates
the number of the step in the cycle over the considered time interval. In each such file, rows
correspond thes cell (on the spacial domain), whereas the columns correspond to the requested for
the output properties, i.e., pH, molality of |H+|, |Ca++|, |Mg++|, |HCO3-|, |CO2(aq)|, as well as
phase volume of calcite and dolomite.

.. _ChemicalState: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html