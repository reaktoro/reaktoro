.. |degC| replace:: °C
.. |H2O| replace:: H\ :sub:`2`\ O
.. |H2Og| replace:: H\ :sub:`2`\ O\ (g)
.. |CO2| replace:: CO\ :sub:`2`
.. |CO2g| replace:: CO\ :sub:`2`\ (g)
.. |MgCl2| replace:: MgCl\ :sub:`2`
.. |CaCl2| replace:: CaCl\ :sub:`2`
.. |%vol| replace:: %\ :sub:`vol`
.. |SiO2| replace:: SiO\ :sub:`2`
.. |CaCO3| replace:: CaCO\ :sub:`3`
.. |NaCl| replace:: NaCl
.. |CaMg(CO3)2| replace:: CaMg(CO\ :sub:`3`)\ :sub:`2`

.. |m2| replace:: m\ :sup:`2`
.. |m3| replace:: m\ :sup:`3`

.. |H+| replace:: H\ :sup:`+`
.. |Ca++| replace:: Ca\ :sup:`2+`
.. |Mg++| replace:: Mg\ :sup:`2+`
.. |HCO3-| replace:: HCO\ :sub:`3`\ :sup:`-`
.. |CO2(aq)| replace:: CO\ :sub:`2` (aq)

.. |slop98| replace:: :download:`slop98.dat <../../../../databases/supcrt/slop98.dat>`

.. |10e-9| replace:: 10\ :sup:`-9`
.. |10e-5| replace:: 10\ :sup:`-5`

Coupling Reaktoro into other reactive transport codes
=====================================================

In this tutorial, we show how Reaktoro can be used in other codes for reactive
transport modeling. Here, you will find that we have split the mass transport
and chemical reaction calculations so that you can clearly see how a dedicated
and advanced transport solver could be combined with Reaktoro's solvers for
chemical reaction calculations. A basic transport solver is used here instead
for the sake of software coupling demonstration.

The reactive transport problem we consider here is also relatively simple,
similar to the one presented in the tutorial
:any:`reactivetransportsolver-calcite-brine`.
Unlike other tutorials, we proceed here first with a step-by-step
explanation of the full script found at the end of this tutorial.

Importing python packages
-------------------------

First, we need to import a few Python packages to enable us to perform the
numerical calculations and plotting.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 1
    :end-before: Step 2

We import the **reaktoro** Python package so that we can use its classes and
methods for performing chemical reaction calculations, **numpy** for working
with arrays, **matplotlib** for plotting capabilities, **joblib** for simple
parallel computing, and, finally, **os**, to provide a portable way of using
operating system dependent functionality.

Initializing auxiliary time-related constants
---------------------------------------------

In this step, we initialize auxiliary time-related constants from seconds up to
years used in the rest of the code.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 2
    :end-before: Step 3

Defining parameters for the reactive transport simulation
---------------------------------------------------------

Next, we define reactive transport and numerical discretization parameters.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 3
    :end-before: Step 4

We specify the considered rock domain by setting coordinates of its left and
right boundaries to 0.0 m and 100.0 m, respectively. The discretization parameters,
i.e., the number of cells and steps in time, are both set to 100. The reactive
transport modeling procedure assumes a constant fluid velocity of 1 m/day
(1.16 · |10e-5| m/s) and the same diffusion coefficient of |10e-9| m2/s for all
fluid species (without dispersivity). The size of the time-step is set to 10
minutes. Temperature and pressure are set to 60 |degC| and 100 bar, respectively,
throughout the whole tutorial.
The chemical equilibrium calculations performed in each mesh cell, at every
time step, using *Gibbs energy minimisation* algorithm (provided by the class
`EquilibriumSolver`_).

The boolean variable ``dirichlet`` is set to ``True`` or ``False`` depending on
which boundary condition is considered in the numerical calculation. Set to
``False`` for imposing the flux of the injected fluid, otherwise, set to
``True`` for imposing the composition of the fluid on the left boundary.

Specifying the quantities and properties to be outputted
--------------------------------------------------------

Before running the reactive transport simulations, we specify the list of
parameters we are interested in outputting. In this case, it is pH, molality of
|H+|, |Ca++|, |Mg++|, |HCO3-|, |CO2(aq)|, as well as a phase volume of calcite
and dolomite.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 4
    :end-before: Step 7

Organizing the main function in three parts
-------------------------------------------

Below is the main function in our script, which shows three parts, each
represented by a Python function documented in the next sections:

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 5

Here, we first create the required folders for the results, we then run the
reactive transport simulation, and, finally, plot the outputted results.

Creating folders for the outputted results
------------------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 6
    :end-before: Step 5

Using **os** package, we create required folders for outputting the obtained
results and for the plot and video files later.

Perform the reactive transport simulation
-----------------------------------------

The reactive transport simulation is performed in the function ``simulate``
(you may want to have a quick look in the code below and move to a
detailed explanation in the subsequent sections):

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7
    :end-before: Step 10

We now explain step-by-step the function ``simulate``. We begin by defining the
chemical system using the `ChemicalEditor`_ class:

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.1
    :end-before: Step 7.2

We specify aqueous and mineral phases that should be considered in the chemical
system. For performance reasons, the aqueous phase is defined by manually
specifying the chemical species, instead of automatic collection of species
from database. There are three pure mineral phase considered: quartz (|SiO2|),
calcite (|CaCO3|), and dolomite (|CaMg(CO3)2|).

Creating the chemical system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step is where we create an object of class `ChemicalSystem`_ using the
chemical system definition details stored in the object ``editor``.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.2
    :end-before: Step 7.3

Initial condition (IC) of the reactive transport problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After constructing the chemical system of interest, we can proceed to the
definition of a chemical equilibrium problem to set the **initial condition**
for the fluid and rock composition that is consistent with our intention of
modeling reactive transport in a porous rock column composed of quartz and
calcite, at 60 |degC| and 100 bar, with a 0.7 NaCl molal brine in equilibrium
with the rock minerals. To ensure equilibrium with both quartz and calcite, the
equilibrium problem setup below considers a relatively large amount for these
minerals (10 moles each).

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.3
    :end-before: Step 7.4

We will later scale the volumes of the aqueous and mineral phases so that they
are consistent with a 10 % porosity and the required volume percentages of the
rock minerals (98 |%vol| of quartz and 2 |%vol| of calcite).

Boundary condition (BC) of the reactive transport problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the boundary condition, we need to specify the composition of the fluid
that is injected into the rock. This is done below, by defining an equilibrium
problem that will later be solved to produce an aqueous fluid with 0.9 molal
NaCl, 0.05 molal |MgCl2|, 0.01 |CaCl2|, and 0.75 molal |CO2|, in a state very
close to |CO2| saturation.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.4
    :end-before: Step 7.5

The composition of the injected aqueous fluid results from the solution of the
above equilibrium problem, in which 1 kg of |H2O| is mixed with 0.90 moles of
|NaCl|, 0.05 moles of |MgCl2|, 0.01 moles of |CaCl2|, and 0.75 moles of |CO2|.

Calculate the equilibrium states for the initial and boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.5
    :end-before: Step 7.6

In this step, we use method ``equilibrate`` to calculate the chemical
equilibrium state of the system with the given initial and boundary equilibrium
conditions stored in the objects ``problem_ic`` and ``problem_bc``. The
numerical solution of each problem results in the objects ``state_ic``
and ``state_bc``, respectively, of class `ChemicalState`_, which stores the
temperature, pressure, and the amounts of every species in the system.


Scaling the volumes of the phases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We scale the volumes of the phases in the initial condition according to the
condition that porosity is 10 % (ratio of fluid volume to total volume), and
the volume fractions of minerals are 98 |%vol| of quartz (|SiO2|) and 2 |%vol|
of calcite (ratio of mineral volume to solid volume).

For the chemical state representing the boundary condition for the injected
fluid composition, we scale its volume to 1 |m3|. This is done so that the
amounts of the species in the fluid is consistent with a mol/|m3| scale.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.6
    :end-before: Step 7.7


Partitioning fluid and solid species
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Only species in fluid phases are mobile and transported by advection and
diffusion mechanisms. The solid phases are immobile.

The code below identify the indices of the fluid and solid species, and also
create arrays to keep track of the amounts of elements in the fluid and solid
partition (i.e., the amounts of elements among all fluid phases, here only an
aqueous phase, and the amounts of elements among all solid phases, here the
mineral phases).

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.7
    :end-before: Step 7.8

We use methods ``indicesFluidSpecies`` and ``indicesSolidSpecies`` of class
`ChemicalSystem`_ to get the indices of the fluid and solid species, which are
then stored in the lists ``ifluid_species`` and ``isolid_species``,
respectively. Then, we define the arrays ``b``, ``bfluid``, ``bsolid``, that
will store, respectively, the concentrations (mol/|m3|) of each element in the
system, in the fluid partition, and in the solid partition at every time step.

The array ``b`` is initialized with the concentrations of the elements at the
initial chemical state, ``state_ic``, using method ``elementAmounts`` of class
`ChemicalState`_. The array ``b_bc`` stores the concentrations of each element
on the boundary.


Create a list of chemical states for the mesh cells
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every mesh cell needs a `ChemicalState`_ object, which is done below with the
creation of a list ``states`` of length ``ncells``, initialized with ``state_ic``
for each cell (i.e., all cells have initially the same fluid and rock composition).

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.8
    :end-before: Step 7.9

Cell coordinates and lengths
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, we generate the coordinates of the mesh nodes (array ``x``) by equally
dividing the interval *[xr, xl]* into ``ncells``. The length between each
consecutive mesh nodes is computed and stored in ``dx`` (the length of the mesh
cells).

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.9
    :end-before: Step 7.10


Creating the equilibrium solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the repeated equilibrium calculation, we define an equilibrium solver
object using `EquilibriumSolver`_  class.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.10
    :end-before: Step 7.11


Define an auxiliary function for creating the output file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The auxiliary function ``outputstate`` uses the class `ChemicalOutput`_ to
output quantities and properties at each mesh cell, at every time step.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.11
    :end-before: Step 7.12

Running the reactive transport simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before proceeding to the reactive transport calculations, we set the initial
time and a counter for the number of time steps.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.12
    :end-before: Step 10

The reactive transport simulation continues until the maximum number of time
steps is achieved. At each time step, we output the progress of the simulation
using function ``outputstate``. Then, we collect the amounts of elements from
fluid and solid partition. We now apply an *operator splitting procedure*, in
which we first update the amounts of elements in the fluid partition
(``bfluid``) using the transport equations (without reactions). These updated
amounts of elements in the fluid partition are now summed with the amounts of
elements in the solid partition (``bsolid``, which remained constant during the
transport step), and thus updating the amounts of elements in the chemical
system (``b``). We now use these updated amounts of elements in the cell
(``b``) to evaluate its new chemical equilibrium state, thus producing new
amounts of the species in both the fluid and solid phases (available in the
list ``states`` of `ChemicalState`_ objects). This chemical reaction
equilibrium calculation step, at each mesh cell, permits, for example,
aqueous species and minerals to react, and thus causes mineral dissolution or
precipitation, depending on how much the amount of mineral species changes. This
can then be used, for example, to compute new porosity value for the cell.


Solving the transport equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reactive transport calculations involve the solution of a system of
advection-diffusion-reaction equations. Due to the use of operator splitting
scheme in this tutorial, we solve the transport equations and then perform the
chemical reaction calculations for each cell before proceeding to the next time
step.

The ``transport`` function below is responsible for solving an
advection-diffusion equation, that is later applied to transport the
concentrations (mol/|m3|) of elements in the fluid partition (*a simplification
that is possible because of common diffusion coefficients and velocities of the
fluid species, otherwise the transport of individual fluid species would be
needed*).

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.10.1
    :end-before: Step 6

The function ``transport`` expects a conservative property (argument ``u``)
(e.g., the concentration (mol/|m3|) of *j*-th element in the fluid given by
``bfluid[j]``), the time step (``dt``), the mesh cell length (``dx``), the
fluid velocity (``v``), the diffusion coefficient (``D``), and the boundary
condition of the conservative property (``g``) (e.g., the concentration of the
*j*-th element in the fluid on the left boundary).

The transport equations are solved with a finite volume method. Its
discretization in space and time (implicit) results in the constants ``alpha``
and ``beta``. These correspond to the diffusion and advection terms in the
equation: ``D*dt/dx**2`` and ``v*dt/dx``, respectively.

Arrays ``a``, ``b``, ``c`` are the diagonals in the tridiagonal matrix that
results by writing all discretized equations in a matrix equation. This system
of linear equations is solved by the tridiagonal matrix algorithm, also known
as the Thomas algorithm:

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 7.10.2
    :end-before: Step 7.10.1


Plotting of the results
-----------------------

The last block of the main routine is plotting of the results.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step 8
    :end-before: Step 7.10.2


The full script for the reactive transport calculation is shown below.

.. literalinclude:: ../../../../demos/python/demo-reactive-transport-calcite-dolomite.py
    :start-after: Step


.. todo::
    Add figures and videos here in this tutorial.

Have you got an issue?
----------------------

Have you found any issue or error in this tutorial? Do you have any
recommendations or you think something is not clear enough? Please, let us know
by filling a new issue here:

.. centered::
    `Reaktoro's GitHub Issues`_

You'll need a GitHub account - but this is easy to sort out if you don't have
one yet!


.. _ChemicalEditor: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalEditor.html
.. _ChemicalSystem: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalSystem.html
.. _ChemicalState: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html
.. _ReactiveTransportSolver: https://reaktoro.org/cpp/classReaktoro_1_1ReactiveTransportSolver.html
.. _ChemicalOutput: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html
.. _KineticPath: https://reaktoro.org/cpp/classReaktoro_1_1KineticPath.html
.. _KineticSolver: https://reaktoro.org/cpp/classReaktoro_1_1KineticSolver.html
.. _EquilibriumSolver: https://reaktoro.org/cpp/classReaktoro_1_1EquilibriumSolver.html
.. _Reaktoro's GitHub Issues: https://github.com/reaktoro/reaktoro/issues/new
