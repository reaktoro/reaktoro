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
.. |CaMg(CO3)2| replace:: CaMg(CO\ :sub:`3`)\ :sub:`2`

.. |H+| replace:: H\ :sup:`+`
.. |Ca++| replace:: Ca\ :sup:`2+`
.. |Mg++| replace:: Mg\ :sup:`2+`
.. |HCO3-| replace:: HC0\ :sub:`3`\ :sup:`-`
.. |CO2(aq)| replace:: C0\ :sub:`2` (aq)

.. |slop98| replace:: :download:`slop98.dat <../../../../databases/supcrt/slop98.dat>`

.. |10e-9| replace:: 10\ :sup:`-9`
.. |10e-5| replace:: 10\ :sup:`-5`

Reactive transport using python tools
=====================================

In this tutorial, we show how Reaktoro can be used in the context of reactive
transport modelling, in which a fluid is continuously injected on the left side
of a rock core. The goal of this tutorial is to illustrate how Reaktoro can be
coupled with other reactive transport simulators, since, unlike the tutorial
:ref:`Reactive transport using  ReactiveTransportSolver`, we provide here the
detailed implementation of the solution of the reactive transport equations,
without using the class `ReactiveTransportSolver`_.

We proceed first with the step-by-step explanation of the script that can be
found in a full length at the very end.

Importing python packages
-------------------------

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 1
    :end-before: Step 2

Here, we import the **reaktoro** Python package so that we can use its classes
and methods for performing the chemical reaction calculations. In addition to
it, we need to import numpy package for working with arrays, **matplotlib** for
plotting capabilities, **joblib** package that enable the set of tools to provide
lightweight pipelining in Python (simple parallel computing), and, finally,
**os** that provides a portable way of using operating system dependent
functionality.

Initializing auxiliary time-related constants
---------------------------------------------

In this step

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 2
    :end-before: Step 3

we initialize auxiliary time-related constants from seconds up to years.

Defining parameters for the reactive transport simulation
---------------------------------------------------------

Next, we define reactive transport and numerical discretization parameters.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 3
    :end-before: Step 4

We specify the considered rock domain by setting coordinates of its left and
right boundaries to 0.0 and 100.0, respectively. The discretization parameters,
i.e., the number of cells and steps in time, are both set to 100. The reactive
transport modelling procedure assumes a constant fluid velocity *v* = 1 *m/day*
(1.16 · |10e-5| *m/s*) and the same diffusion coefficient *D* = |10e-9| *m2/s*
for all fluid species without dispersivity. The size of the time-step is set to
10 minutes. The temperature is set to 60 |degC|, whereas the pressure is fixed
to 100 bar. Parameter ``dirichlet`` is initialized by **True** or **False**
depending on which boundary conditions are considered in the numerical test. If
it is set to **False**,  the natural flux conditions are imposed, otherwise, the
Dirichlet one is selected. The last parameter ``smrt_solv`` is initialized by
**True** or **False** depending on which equilibrium solver is used. Setting it
to **True** enables the usage of the `SmartEquilibriumSolver`_, which exploits
*smart on-demand learning* of entirely new states of the chemical system.

Defining the list of the output quantities
------------------------------------------

Before running the reactive transport simulations, we specify the list of
parameters we are interested in outputting. In this case, it is pH, molality of
|H+|, |Ca++|, |Mg++|, |HCO3-|, |CO2(aq)|, as well as a phase volume of calcite
and dolomite.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 4
    :end-before: Step 7

Specifying the structure of the numerical experiment
----------------------------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 5

To provide the entry point of the python script, we use
``if __name__ == '__main__':``. The structure of the numerical experiment contains
three main blocks, i.e., creating the required folders for the results, running
the simulations, and plotting the obtained results.

Creating folders for the results' output
----------------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 6
    :end-before: Step 5

Using **os** package, we create required folders for outputting the obtained
results.

Perform the reactive transport simulation
-----------------------------------------

The reactive transport simulation is presented by the corresponding method
``simulate()``.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7
    :end-before: Step 10

Its content will be explained below step-by-step. We begin by defining the
chemical system using the`ChemicalEditor`_ class:

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.1
    :end-before: Step 7.2

We specify **aqueous** and **mineral** phases that should be considered in the
chemical system. For performance reasons, the aqueous phase is defined by
manually specifying the chemical species, instead of automatic collection of
the collection of species from database. There are three pure mineral phase
considered: quartz (|SiO2|), calcite (|CaCO3|), and dolomite (|CaMg(CO3)2|).

Creating the chemical system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step is where we create an object of class `ChemicalSystem`_ using the
chemical system definition details stored in the object ``editor``.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.2
    :end-before: Step 7.3

Initial condition (IC) of the reactive transport problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After constructing the chemical system of interest, we can proceed to the
*definition of the chemical reaction problem*. First, we specify its **initial
condition** with given conditions for *temperature*, *pressure*, and *amounts
of elements* that are consistent with an intention of modelling reactive transport
of injected |NaCl|-|MgCl2|-|CaCl2| brine into the rock-fluid composition of
quartz and calcite at 60 |degC| and 100 bar. In particular, we consider a 0.7
molal |NaCl| brine in equilibrium with the rock minerals (with a calculated pH
of 10.0).

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.3
    :end-before: Step 7.4

Boundary condition (BC) of the reactive transport problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we define the **boundary condition** of the constructed chemical system
with its *temperature*, *pressure*, and *amounts of elements*.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.4
    :end-before: Step 7.5

In particular, we prescribe the amount of injected aqueous fluid resulting from
the mixture of 1 kg of water with 0.90 moles of |NaCl|, 0.05 moles of |MgCl2|,
0.01 moles of |CaCl2|, and 0.75 moles of |CO2|, in a state very close to |CO2|
saturation. The temperature and the pressure stay the same, i.e., 60 |degC| and
100 bar, respectively.

Calculate the equilibrium states for the initial and boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.5
    :end-before: Step 7.6

In this step, we use ``equilibrate()`` method to calculate the chemical
equilibrium state of the system with the given initial and boundary equilibrium
conditions stored in the object ``problem_ic`` and ``problem_bc``. For this
calculation, an efficient **Gibbs energy minimization** computation is used to
determine the species' amounts that correspond to a state of minimum Gibbs energy
in the system, while satisfying the prescribed conditions for temperature,
pressure, and element amounts. The result is stored in the objects ``state_ic``
and ``state_bc`` of class `ChemicalState`_.

.. attention::

    In the vast majority of cases, there is one object of class `ChemicalSystem`_
    in the code and one or more objects of class `ChemicalState`_ describing
    different states of the chemical system. Reaktoro differentiates these two
    independent concepts: *chemical system definition* and *chemical system state*.


Scaling the volumes of the phases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We scale the volumes of the phases in the initial condition according to the
following composition, i.e., 97.73 |%vol| |SiO2| (quartz) and 2.27 |%vol| |CaCO3|
(calcite) with the porosity of 10%. Moreover, the volume of the fluid in the
boundary condition state is scaled to 1 m3. This is done so that the amounts of
the species in the fluid is consistent with a mol/m3 scale.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.6
    :end-before: Step 7.7

Specifying discretization structures needed for the reactive transport
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the next block of code, we provide structures needed for the discretization
of the problem and reactive transport simulation.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.7
    :end-before: Step 7.8

We start by fetching the indices of the fluid and solid species to the
``ifluid_species`` and ``isolid_species`` lists, respectively. Then, we define
the numpy arrays ``b``, ``bprev``, ``bfluid``, ``bsolid``, that will store the
concentrations of each element in each mesh cell for the whole system at the
current time step as well as its fluid and solid partitions.

The array corresponding to the current step of the reactive transport is
initialized by the concentrations of elements at the initial chemical state using
``state_ic.elementAmounts()``. The array ``b_bc`` will store the concentrations
of each element on the boundary. The list ``states`` (of length ``ncells``)
contains chemical states for each cell.

Next, we generate the set of degrees of freedom (DOFs) discretizing the domain
represented by the interval *[xr, xl]* (can be regarded as a mesh). The
corresponding spacial mesh size is defined by ``dx``.


Creating the equilibrium solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the repeated equilibrium calculation, we define an equilibrium solver object
using either `EquilibriumSolver`_ or `SmartEquilibriumSolver`_ classes, which
are initialized by a considered chemical system. Here, if the parameter
``smrt_solv`` is set to **False**, the classical approach (performed by the class
`EquilibriumSolver`_) is considered. Otherwise, class `SmartEquilibriumSolver`_
is used.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.8
    :end-before: Step 7.9

In reactive transport simulations, we have to perform many equilibrium
calculations in sequence. The usage of `EquilibriumSolver`_ is more appropriate
because each call of ``equilibrate()`` method has a computational overhead of
creating a new object of class `EquilibriumSolver`_. Therefore, we create this
solver object only once and then use it subsequently for all other equilibrium
calculations. Analogously, for repeated chemical kinetics calculations (e.g.,
in a reactive transport simulation, where kinetics is performed for each mesh
cell/node), consider using the class `KineticSolver`_  for avoiding some overhead
of class `KineticPath`_.

Define an auxiliary function for creating the output file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The auxiliary routine ``outputstate()`` is a function that uses the class
`ChemicalOutput`_ to output requested by the user information into the output
file.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.9
    :end-before: Step 7.10

Running the reactive transport simulation loop
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before proceeding to the simulation of reactive transport on the considered
interval, we set the initial time and a counter for the step in this cycle. In
the variable ``ndigits`` we store the number of digits in the selected amount of
steps needed for the naming the files with chemical properties.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.10
    :end-before: Step 10

The cycle for the reactive transport simulation proceeds until we haven't made
all the steps in time. At each time step, we output the progress of the
simulations by auxiliary routine ``outputstate()``. Then, we collect the amounts
of elements from fluid and solid partition. This is done according to the
**operator splitting procedure**. First, we update the amounts of elements in
the fluid partition using the transport equations.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.10.1
    :end-before: Step 6

The transport is performed by the function ``transport`` that expects the amounts
of the *j*th element, current time and space discretization step, parameters of
the reactive transport velocity and diffusion coefficient, and the amounts of
*j*-th element on the boundary. This transport solver method has the following
structure. It defines the discretization constants ``alpha`` and ``beta`` that
correspond to the diffusion and advection terms in the equation, i.e.,
``D*dt/dx**2`` and ``v*dt/dx``, respectively. Arrays ``a``, ``b``, ``c`` follow
from the the finite difference discretization of of the reaction-advection
equation. Finally, the obtained system of equation is solved by the tridiagonal
matrix algorithm, also known as the Thomas algorithm.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 7.10.2
    :end-before: Step 7.10.1

-------------------------------------------------------------------------------
The second step of **operator splitting procedure** is to update the total
amount of elements by accounting the earlier reconstructed fluid partition and
the solid partition that stays constant with respect to time.

 .. code-block:: python

    # Update the amounts of elements in both fluid and solid partitions
        b[:] = bsolid + bfluid

Plotting of the results
-----------------------

The last block of the main routine is plotting of the results.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step 8
    :end-before: Step 7.10.2


The above step-by-step explanation is summarized in a script below.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-at: Step

.. _ChemicalEditor: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html
.. _ChemicalSystem: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html
.. _ChemicalState: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html
.. _ReactiveTransportSolver: https://reaktoro.org/cpp/classReaktoro_1_1ReactiveTransportSolver.html
.. _SmartEquilibriumSolver: https://reaktoro.org/cpp/classReaktoro_1_1SmartEquilibriumSolver.html
.. _ChemicalOutput: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html
.. _KineticPath: https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html
.. _KineticSolver: https://reaktoro.org/cpp/classReaktoro_1_1KineticSolver.html
.. _EquilibriumSolver: https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html