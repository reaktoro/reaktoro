Chemical Equilibrium -- Basic
=============================

We present below a Python script that performs a single-phase chemical
equilibrium calculation using Reaktoro. The single phase in the chemical system
is an *aqueous phase* with its chemical species automatically collected from a
given thermodynamic database. We are interested in calculating the amounts of
every aqueous species after we let 1 kg of H2O, 0.5 mol of NaCl, and 1 mol of
CO2 to react to a state of chemical equilibrium under constant temperature and
pressure of 60 °C and 100 bar respectively.

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py

We explain next every step in this chemical reaction modeling exercise.

Importing the reaktoro Python package
-------------------------------------

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py
    :start-at: Step 1
    :end-before: Step 2

Using Reaktoro in Python requires first an import of the python package
**reaktoro**. From this point on, we are able to use the library components of
Reaktoro (classes, methods, constants), which are needed to define our chemical
system and chemical reaction modeling problems.

.. note::

    To simplify the tutorials, we use ``from reaktoro import *``, which imports
    all components of the **reaktoro** package into the default Python
    namespace, which can potentially create name conflicts. For your
    applications, consider instead using ``import reaktoro as rkt``,
    and then refer to Reaktoro's classes and methods as ``rkt.Database``,
    ``rkt.ChemicalSystem``, ``rkt.equilibrate``, and so forth.

Initializing a thermodynamic database
-------------------------------------

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py
    :start-at: Step 2
    :end-before: Step 3

.. |supcrt98| replace:: :download:`supcrt98.xml <../../databases/supcrt/supcrt98.xml>`
.. |supcrt07| replace:: :download:`supcrt07.xml <../../databases/supcrt/supcrt07.xml>`
.. |supcrt98-organics| replace:: :download:`supcrt98-organics.xml <../../databases/supcrt/supcrt98-organics.xml>`
.. |supcrt07-organics| replace:: :download:`supcrt07-organics.xml <../../databases/supcrt/supcrt07-organics.xml>`
.. |slop98| replace:: :download:`slop98.dat <../../databases/supcrt/slop98.dat>`
.. |slop07| replace:: :download:`slop07.dat <../../databases/supcrt/slop07.dat>`

In this step we initialize a ``Database`` object with the |supcrt98| database
file. This database was generated from the original SUPCRT92 database file
|slop98|. You are welcome to inspect these files and learn more about the
chemical species available in them. You can also read more about the available
thermodynamic databases supported in Reaktoro at :ref:`Thermodynamic
Databases`.

Defining your chemical system
-----------------------------

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py
    :start-at: Step 3
    :end-before: Step 4

Reaktoro is a general-purpose chemical solver that avoids as much as possible
infering specific assumptions about your problems. Thus, you need to specify
how your chemical system should be defined. This encompass the specification of
all phases in the system as well as the chemical species that compose each
phase. By using the ``ChemicalEditor`` class, you can achieve this.

In this step, we create an object of class ``ChemicalEditor`` and specify that
an aqueous phase should be considered in the chemical system with the method
call ``addAqueousPhase``. The string ``'H O Na Cl C'`` tells Reaktoro to go
through the thermodynamic database and collect all aqueous species that can be
created by combining the chemical elements H, O, Na, Cl, and C.

.. note::

    An automatic search for chemical species can result in a large number of
    species in the phase, and thus causing the chemical reaction calculations
    to be more computationally expensive. If you are using Reaktoro in
    demanding applications (e.g., as a chemical solver in a reactive transport
    simulator), you might want to manually specify the chemical species of each
    phase in your chemical system. This can be achieved by providing a list of
    species names as shown below:

    .. code-block:: python

        editor.addAqueousPhase([
            'H2O(l)',
            'H+',
            'OH-',
            'Na+',
            'Cl-',
            'HCO3-',
            'CO3--',
            'CO2(aq)'
            ])

.. caution::

    If you manually specify the chemical species in a phase, you need to make
    sure that they exist in the thermodynamic database with **exact same
    name**! Replacing ``'CO2(aq)'`` above with ``'CO2'`` will cause an error if
    the database has no aqueous species with such name.

Constructing your chemical system
---------------------------------

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py
    :start-at: Step 4
    :end-before: Step 5

This step is where we create an object of class ``ChemicalSystem`` using the
chemical system definition details stored in the object ``editor``.

.. note::

    The ``ChemicalSystem`` is perhaps the main class in Reaktoro. An object of
    this class stores the phases, species, elements, and reactions in the
    chemical system, as wel as provides the means to compute many types of
    thermodynamic properties, such as *standard thermodynamic properties*
    (e.g., standard Gibbs energies, standard enthalpies, and standard volumes
    of species), and *thermo-chemical properties* (e.g., activity and activity
    coefficients of species; density, enthalpy and internal energy of phases).

As you learn more about other Reaktoro's classes, you will note that an object
of class ``ChemicalSystem`` is almost always needed for their initialization!
This is because no chemical reaction calculation can be performed without the
details of the chemical system and the methods to evaluate thermodynamic
properties of its species, phases, and reactions.

Defining the chemical equilibrium problem
-----------------------------------------

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py
    :start-at: Step 5
    :end-before: Step 6

We have now defined and constructed our chemical system of interest, enabling
us to move on to the next step in Reaktoro's modeling workflow: *formulating
and solving interesting chemical reaction problems*. In this tutorial we are
interested in computing the chemical equilibrium state of an aqueous phase for
given equilibrium conditions of:

* *temperature*;
* *pressure*; and
* *amounts of chemical elements*.

You might be wondering if this is correct, because  what we actually provided
above were the amounts of substances H2O, NaCl, and CO2 using the method
``EquilibriumProblem.add``. Yes, but here is what happens behind the scenes:
Reaktoro parses these chemical formula strings to determine the chemical
elements and their coefficients in the formula. Once this is done, the amount
of the chemical elements hidden inside ``EquilibriumProblem`` are incremented
accordingly. For this, the coefficients of the elements in the formula are used
to obtain their amounts from the given amount of the substance.

Thus, the following method call:

.. code-block:: python

    problem.add('CO2', 1.0, 'mol')

is equivalent to:

.. code-block:: python

    problem.add('C', 1.0, 'mol')
    problem.add('O', 2.0, 'mol')

The table below shows the amounts of chemical elements resulting from the
combination of 1 kg of H2O, 1 mol of CO2, and 0.5 mol of NaCl (assume 55.5 mol
of H2O is roughly 1 kg of H2O).

.. table:: Amounts of chemical elements obtained from the
    recipe 1.0 kg of H2O (55.5 mol), 1.0 mol of CO2, and 0.5 mol of NaCl.
    :widths: 1 1
    :align: center

    ======= ============
    Element Amount (mol)
    ======= ============
    H       111.0
    O       57.5
    Na      0.5
    Cl      0.5
    C       1.0
    ======= ============

.. danger::
    Now that you know that an equivalent chemical equilibrium problem could be defined with:

    .. code-block:: python

        problem.add('H' 111.0, 'mol')
        problem.add('O' 57.5, 'mol')
        problem.add('Na' 0.5, 'mol')
        problem.add('Cl' 0.5, 'mol')
        problem.add('C' 1.0, 'mol')

    you might want to *adventure* in manually specifying different values for
    the amounts of elements. Just be extra careful with the values you provide
    as this could result in an *infeasible* chemical equilibrium state. Here is
    a simple example of conditions that result in an infeasible equilibrium
    state. Consider a chemical system containing only a gaseous phase with
    gases H2O(g) and CO2(g). Find non-negative amounts for H2O(g) and CO2(g)
    when the given amounts of chemical elements are: 2 mol of H, 1 mol of O,
    and 1 mol of C.

    Please note that we are not condemning the input form shown above in terms
    of element amounts, but only telling you to be attentive with the values
    you input. If you are using Reaktoro as a chemical reaction solver in a
    reactive transport simulator, for example, you'll most likely need to work
    directly with given amounts of elements, which shows that this input form
    is required in certain cases. For such time-dependent modeling problems, you
    often only need to ensure that the initial conditions for elements amounts
    result in feasible initial species amounts.

    .. todo::
        Give a reference here to a reactive transport tutorial page where this
        is better explained.

.. tip::

    The substance formulas given in the ``EquilibriumProblem.add`` method
    **can, but do not need to**, correspond to names of chemical species in the
    thermodynamic database. Even unusual, if not strange, substance formulas,
    such as HCl3(NaO)4C13, would be understood by that method. *We do not
    promise, however, that you will obtain a feasible chemical equilibrium
    state with unrealistic conditions!*

Calculating the chemical equilibrium state
------------------------------------------

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py
    :start-at: Step 6
    :end-before: Step 7

In this step we use the ``equilibrate`` function to calculate the chemical
equilibrium state of the system with the given equilibrium conditions stored in
the object ``problem``. For this calculation, Reaktoro uses an efficient
**Gibbs energy minimization** computation to determine the species amounts that
correspond to a state of minimum Gibbs energy in the system, while satisfying
the prescribed amount conditions for temperature, pressure, and element
amounts. The result is stored in the object ``state``, of class
``ChemicalState``.

.. attention::

    In the vast majority of cases, you'll only have one object of class
    ``ChemicalSystem`` in your code and one or more objects of
    ``ChemicalState`` describing different states of your chemical system!
    Reaktoro differentiates these two independent concepts: *chemical system
    definition* and *chemical system state*.

.. note::

    The method ``equilibrate`` is a convenient function in Reaktoro. Consider
    using the class ``EquilibriumSolver`` for more advanced requirements. For
    example, if you have to perform many equilibrium calculations in sequence.
    The ``equilibrate`` method has a computational overhead because every call
    creates a new object of class ``EquilibriumSolver``. Preferably, this
    object should be created only once, and then used subsequently for all
    other equilibrium calculations.

Outputting the calculated chemical state to a file
--------------------------------------------------

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py
    :start-at: Step 7
    :end-before: Step 8

You've performed your chemical equilibrium calculation and you now want to
inspect the computed compositional state and its thermodynamic properties. Here
is a convenient wy of doing so: outputting the chemical state to a file, here
named ``result.txt``.

.. role:: python(code)
    :language: python3

.. tip::
    You can also print the chemical state directly in the terminal with:

    .. code-block:: python

        print(state)

    Or in C++:

    .. code-block:: c++

        std::cout << state << std::endl;

Printing the amounts of some aqueous species
--------------------------------------------

.. literalinclude:: ../../demos/python/demo-equilibrium-basic.py
    :start-at: Step 8

Here is just a small demonstration of getting species amount information from a
``ChemicalState`` object using the method ``ChemicalState.speciesAmount`` to
extract the amount of a few chemical species. There are many more methods you
could use, and you are welcome to inspect the interface of this class.

.. todo:: Add reference here to the documentation of ``ChemicalState`` class.
