.. |degC| replace:: °C
.. |m2| replace:: m\ :sup:`2`\
.. |m3| replace:: m\ :sup:`3`\
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

.. |H+| replace:: H\ :sup:`+`
.. |Ca++| replace:: Ca\ :sup:`2+`
.. |Mg++| replace:: Mg\ :sup:`2+`
.. |HCO3-| replace:: HCO\ :sub:`3`\ :sup:`-`
.. |CO2(aq)| replace:: CO\ :sub:`2` (aq)

.. |10e-21| replace:: 10\ :sup:`-21`
.. |10e-9| replace:: 10\ :sup:`-9`
.. |10e-5| replace:: 10\ :sup:`-5`

.. |slop98| replace:: :download:`slop98.dat <../../../../databases/supcrt/slop98.dat>`

Reactive transport of |CO2|-saturated brine along a porous rock column
======================================================================

In this tutorial, we show how Reaktoro can be used for one-dimensional reactive
transport calculations for modeling the geochemical reactions that occur along
a porous rock column as an aqueous fluid is continuously injected on its left
side.

The injected fluid is a brine with 0.9 molal NaCl, 0.05 molal |MgCl2|, 0.01
molal |CaCl2| and almost |CO2|-saturated, with 0.75 molal of |CO2| dissolved.

The porous rock is initially composed of minerals quartz (|SiO2|) and calcite
(|CaCO3|). The initial porosity is 10 %, and the initial volume percentages of
the minerals are: 98 |%vol| of quartz, 2 |%vol| calcite. The initial conditions
for the fluid in the rock is a 0.7 molal |NaCl| brine in equilibrium with the
existing rock minerals calcite and quartz. These initial fluid and rock
composition conditions are uniform throughout the rock core. We assume a rock
column length of 100 m at temperature 60 |degC| and 100 bar throughout.

.. admonition:: Assumptions

    To simplify this tutorial, the following assumptions are made:

    * Chemical equilibrium is used for modeling the chemical reactions
      in this problem, not only for reactions between aqueous-aqueous species,
      but also for those between mineral and aqueous species.
    * A uniform and constant velocity field is imposed and it is not updated
      by solving, for example, Darcy equation.
    * Both temperature and pressure are also kept constant along the rock.


.. todo::

    Create other more complicated tutorials that do not consider the previous
    assumptions so that chemical kinetics is considered for mineral dissolution
    and precipitation reactions and the partial chemical equilibrium assumption
    is assumed for aqueous and gaseous species. Also create tutorials in which
    velocity and pressure fields are physically consistent with each other and
    computed via the solution of Darcy equation for porous media flow coupled
    with mass transport equations.

More details about the problem setup can be found in the script below, and also
in the step-by-step explanation that follows it.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step


Importing the reaktoro package
------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 1
    :end-before: Step 2

First, we import the **reaktoro** Python package so that we can use its classes
and methods for performing the chemical reaction calculations.

Defining auxiliary time-related constants
-----------------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 2
    :end-before: Step 3

In this step we initialize auxiliary time-related constants from seconds to
years. This is only done for convenience, so that we can specify later, for
example, fluid velocity as 1 m/day.

Defining parameters for the reactive transport simulation
---------------------------------------------------------

Next, we define reactive transport and numerical discretization parameters.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 3
    :end-before: Step 4

First, we define the considered rock domain by setting coordinates of its left
and right boundaries to 0.0 m and 100.0 m, respectively. The number of cells
and steps in time are both set to 100. The reactive transport modelling
procedure assumes a constant fluid velocity of *v* = 1 m/day (1.16 · |10e-5|
m/s) and the same diffusion coefficient *D* = |10e-9| |m2|/s for all fluid
species. The size of the time-step is set to be half of day. The temperature
and pressure of the fluids are selected to be 60 |degC| and 100 bar,
respectively.

Defining the chemical system
----------------------------

We need to define a chemical system that can represent both our fluid and rock.
We use class `ChemicalEditor`_ below to define a system with an aqueous phase
and three mineral phases: quartz, calcite and dolomite. Initially, our rock has
no dolomite (|CaMg(CO3)2|), but since this is a mineral that could potentially
precipitate given the fluid composition injected (containing |CaCl2| and
|MgCl2| dissolved), we add it here in the chemical system to ensure that the
calculations are able to model dolomite precipitation.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 4
    :end-before: Step 5

.. note::
    The aqueous phase is defined above by using a list of compounds, which is
    then broken automatically by Reaktoro into a list of element names. These
    element names are then used to find in the database all the aqueous species
    that could be formed out of them.

Constructing the chemical system
--------------------------------

This step is where we create an object of class `ChemicalSystem`_ using the
chemical system definition details stored in the object ``editor``.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 5
    :end-before: Step 6


Initial condition (IC) for the fluid composition
---------------------------------------------------------------

Below, we define the **initial condition** for the fluid composition in the
rock. We want an aqueous fluid that is 0.7 molal of NaCl and in equilibrium
(saturated) with calcite (|CaCO3|) and quartz (|SiO2|). To achieve this, we mix
1 kg of |H2O|, 0.7 mol of NaCl, and plenty of calcite and quartz (10 mol each) to
ensure that the aqueous solution is saturated with respect to these minerals.
Temperature and pressure are set to 60 °C and 100 bar respectively.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 6
    :end-before: Step 7


Boundary condition (BC) for the fluid composition
-------------------------------------------------

Next, we define the **boundary condition** for the fluid composition on the
left side of the rock, which should be the one that represents the fluid being
continuously injected: 0.9 molal NaCl, 0.05 molal |MgCl2|, 0.01 molal |CaCl2|
and almost |CO2|-saturated, with 0.75 molal of |CO2| dissolved. To achieve
this, we mix 1 kg of |H2O| with 0.9 mol of NaCl, 0.05 mol of |MgCl2|, 0.01 mol
of |CaCl2|, and 0.75 mol of |CO2|. Temperature and pressure are also set to 60
°C and 100 bar respectively.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 7
    :end-before: Step 8


Calculating the IC and BC fluid compositions
--------------------------------------------

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 8
    :end-before: Step 9

In this step, we use the ``equilibrate`` function to calculate the chemical
equilibrium state of the system with the given initial and boundary equilibrium
conditions stored in the object ``problem_ic`` and ``problem_bc``.  The result
is stored in the objects ``state_ic`` and ``state_bc`` of class
`ChemicalState`_.


Scaling the phases in the initial condition
-------------------------------------------

The initial chemical state ``state_ic`` computed before has, at this point,
phases with volumes that do not correspond to our desired porosity of 10 % and
rock mineral composition of 98 |%vol| of quartz and 2 |%vol| of calcite.

To obtain this, we scale the volumes of the aqueous and mineral phases as shown
below:

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 9
    :end-before: Step 10

.. note::

    After this scaling step, the sum of the phase volumes in ``state_ic`` is 1
    |m3|. This also ensures that the amounts of the species in the chemical
    system are normalized by |m3|, and thus they can be regarded as
    concentrations in unit of mol/|m3| (*bulk volume, not fluid volume!*).


Scaling the boundary condition state
------------------------------------

Next, we scale the boundary condition state to 1 |m3|, so that we have the
amounts of fluid species in ``state_bc`` also normalized by |m3|.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 10
    :end-before: Step 11

.. note::
    The chemical state represented by ``state_bc`` has no other stable phase
    than the aqueous phase (i.e., all mineral phases have zero or negligible
    amounts such as |10e-21| mol).

Creating the mesh
-----------------

We define the mesh with the class `Mesh`_ in order to use in the initialization
of class `ReactiveTransportSolver`_ later.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 11
    :end-before: Step 12

Here, we specify the number of cells in the mesh and the x-coordinates of the
left and right boundaries (in m).

Creating a chemical field object
--------------------------------

We have been using class `ChemicalState`_ to represent an individual chemical
state. We will now use class `ChemicalField`_ to represent a collection of
chemical states: one for each mesh cell.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 12
    :end-before: Step 13

.. note::
    Different initial conditions across the mesh cells is possible by assigning
    different chemical states to each mesh cell. Here, the same chemical state
    in ``state_ic`` is used for all cells.


Initializing the reactive transport solver
------------------------------------------

At last, we define the object responsible for the solving the reactive
transport problem, which is handled by the class `ReactiveTransportSolver`_:.:

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 13
    :end-before: Step 14

Here, we set the mesh and problem parameters such as velocity, diffusion
coefficient, the chemical state representing the boundary condition, and the
time step. We also initialize the reactive solver object with the chemical
field object specified on the previous step, at this point containing the
initial condition for the chemical state of each mesh cell.

Defining the output quantities
------------------------------

Before starting the reactive transport calculations, we define the quantities
that will be output for every mesh cell, at every time step. For this, we use
an object of the class `ChemicalOutput`_:

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 14
    :end-before: Step 15

The name of the output file is to ``reactive-transport.txt``. We specify the
parameters that we are interested in outputting. In this case, it is pH,
molality of |H+|, |Ca++|, |Mg++|, |HCO3-|, |CO2(aq)|, as well as a phase volume
of calcite and dolomite.


Running the reactive transport simulation
-----------------------------------------

As shown below, we perform a sequence of reactive transport calculations, one
for each time step, during which the chemical state of each mesh cell is
updated. The iterations continue until the maximum number of steps is
achieved.

.. literalinclude:: ../../../../demos/python/demo-reactivetransportsolver-calcite-brine.py
    :start-at: Step 15

At each time step, we print the progress of the simulation. We then use method
``step`` of class `ReactiveTransportSolver`_ to perform a single reactive
transport time-stepping. This method also produces a new output file containing
the requested output properties for every mesh cell. In each such file, rows
correspond to cells, whereas the columns correspond to the requested output
properties, i.e., pH, molality of |H+|, |Ca++|, |Mg++|, |HCO3-|, |CO2(aq)|, as
well as the phase volume of calcite and dolomite.

.. todo::

    Use the result files to general plots and videos, and add them here in the
    tutorial.


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
.. _TransportSolver: https://reaktoro.org/cpp/classReaktoro_1_1TransportSolver.html
.. _SmartEquilibriumSolver: https://reaktoro.org/cpp/classReaktoro_1_1SmartEquilibriumSolver.html
.. _ChemicalOutput: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalOutput.html
.. _Mesh: https://reaktoro.org/cpp/classReaktoro_1_1Mesh.html
.. _ChemicalField: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalField.html
.. _Reaktoro's GitHub Issues: https://github.com/reaktoro/reaktoro/issues/new
