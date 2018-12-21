.. |H2O| replace:: H\ :sub:`2`\ O
.. |H2Og| replace:: H\ :sub:`2`\ O\ (g)
.. |CO2| replace:: CO\ :sub:`2`
.. |CO2g| replace:: CO\ :sub:`2`\ (g)
.. |MgCl2| replace:: MgCl\ :sub:`2`
.. |CaCl2| replace:: CaCl\ :sub:`2`

Solubility of |CO2| in NaCl brines
==================================

In this tutorial, we show how Reaktoro can be used to compute the solubility of
|CO2| in a 1 molal NaCl brine at temperature 60 °C and pressure 100 bar. We
show no *magical function* to perform such calculation in a single line of
code, but instead a sequence of steps using Reaktoro's components (classes,
methods) to enrich your understanding of how Reaktoro can be used for solving
this and many other different chemical reaction modeling problems.

To calculate the solubility of |CO2| in a NaCl brine, we need two phases in our
chemical system: an *aqueous phase* to represent our NaCl brine, and a *gaseous
phase* to represent our |CO2| gas. Next, we formulate and solve a *chemical
equilibrium problem*, in which 10 mol of |CO2| is mixed with 1 kg of |H2O| and
1 mol of NaCl, at 60°C and 100 bar. The solution to this problem is a *chemical
equilibrium state* from which we can inspect how much |CO2| exist in the
gaseous phase and compare this with the initial amount of |CO2| we used to mix
with |H2O| and NaCl.

.. note::

    Our **1 molal NaCl brine** is here represented by the mixture of 1 kg of |H2O|
    and 1 mol of NaCl.

.. note::

    Given that the mole mass of |CO2| is roughly 44 g/mol, 10 mol of |CO2| is
    approximately 440 g! Mixing this amount of |CO2| with 1 kg of |H2O| will most
    likely result in a chemical equilibrium state in which the aqueous phase is
    saturated with |CO2| and the remaining |CO2| exists as a gas (or super-critical
    fluid depending on the temperature and pressure). If we choose a small
    amount for |CO2|, we can end up with an equilibrium state in which all |CO2|
    has been dissolved in the aqueous phase, and we will then not be able to
    determine its solubility.

We present below the Python script that performs a multi-phase, multi-species
chemical equilibrium calculation using Reaktoro to determine the solubility of
|CO2| in a 1 molal NaCl brine at temperature 60 °C and pressure 100 bar.

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step

You find next a step-by-step explanation of the above script.

Importing the reaktoro Python package
-------------------------------------

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
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

Thermodynamic databases are essential for modeling chemically reactive systems
using Reaktoro. We need a database from where we collect data of substances
that will compose our phases of interest in a multi-phase chemical system. A
thermodynamic database also contains model parameters for the evaluation of
standard thermodynamic properties of species and/or reactions (e.g., standard
Gibbs energies, equilibrium constants).

In this step:

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step 2
    :end-before: Step 3

.. |supcrt98| replace:: :download:`supcrt98.xml <../../../../databases/supcrt/supcrt98.xml>`
.. |slop98| replace:: :download:`slop98.dat <../../../../databases/supcrt/slop98.dat>`

we initialize a ``Database`` object with the |supcrt98| database
file. This database was generated from the original SUPCRT92 database file
|slop98|. You are welcome to inspect these files and learn more about the
chemical species available in them. You can also read more about the available
thermodynamic databases supported in Reaktoro at :ref:`Thermodynamic
Databases`.

.. note::

    All SUPCRT92 thermodynamic databases have been embedded into Reaktoro.
    Thus, you don't actually need to have a database file named
    ``supcrt98.xml`` in a local directory when initializing the ``Database``
    object. If you want to use a customized database file, however, also named
    ``supcrt98.xml``, then your local file will be used instead.

.. tip::

    If you are using a customized version of a thermodynamic database, consider
    changing its name (e.g., ``custom-supcrt98.xml``) to avoid an accidental
    use of an embedded database. This can happen if you do not give a correct
    path to your custom database file.


Defining the chemical system
----------------------------

Reaktoro is a general-purpose chemical solver that avoids as much as possible
presuming specific assumptions about your problems. Thus, you need to specify
how your chemical system should be defined. This encompass the specification of
all phases in the system as well as the chemical species that compose each
phase. By using the ``ChemicalEditor`` class, you can conveniently achieve
this as shown below:

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step 3
    :end-before: Step 4

In this step, we create an object of class ``ChemicalEditor`` and specify that
an **aqueous phase** and a **gaseous phase** should be considered in the
chemical system. The aqueous phase is defined such that its species are all
those aqueous species in the selected thermodynamic database that can be
created by combining the chemical elements H, O, Na, Cl, and C. The gaseous
phase is defined with only one gaseous species: |CO2g|.

.. note::

    The selected elements H, O, Na, Cl, and C in the definition of the aqueous
    phase and the choice of |CO2g| as the single gaseous species composing the
    gaseous phase consistently represent our intentions of calculating the
    solubility of |CO2| in a NaCl brine. If we decide to do a similar
    computation, but with a NaCl-|MgCl2|-|CaCl2| brine, then that original list
    of elements would need to be incremented with Mg and Ca.


.. note::

    An automatic search for chemical species can result in a large number of
    species in the phase, potentially causing the chemical reaction
    calculations to be more computationally expensive. If you are using
    Reaktoro in demanding applications (e.g., as a chemical solver in a
    reactive transport simulator), you might want to manually specify the
    chemical species of each phase in your chemical system. This can be
    achieved by providing a list of species names as shown below:

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

    This is exactly what we did for the definition of the gaseous phase. If we
    had done instead:

    .. code-block:: python

        editor.addGaseousPhase('C O')

    then other gases would be considered, such as CO(g) and O2(g), which are
    not of interest in our modeling problem.

.. caution::

    If you manually specify the chemical species in a phase, you need to make
    sure that they exist in the thermodynamic database with **exact same
    name**! Replacing ``'CO2(g)'`` above with ``'CO2'`` will cause an error if
    the database has no gaseous species with such name.

Constructing the chemical system
--------------------------------

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step 4
    :end-before: Step 5

This step is where we create an object of class ``ChemicalSystem`` using the
chemical system definition details stored in the object ``editor``.

.. note::

    ``ChemicalSystem`` is perhaps the main class in Reaktoro. An object of this
    class stores the phases, species and elements in our defined chemical
    system, as wel as provides the means to compute many types of thermodynamic
    properties, such as *standard thermodynamic properties* (e.g., standard
    Gibbs energies, standard enthalpies, and standard volumes of species), and
    *thermo-chemical properties* (e.g., activity and activity coefficients of
    species; density, enthalpy and internal energy of phases). As you learn
    more about other Reaktoro's classes, you will note that an object of class
    ``ChemicalSystem`` is almost always needed for their initialization!

Defining the chemical equilibrium problem
-----------------------------------------

We have now defined and constructed our chemical system of interest, enabling
us to move on to the next step in Reaktoro's modeling workflow: *defining our
chemical reaction problems*. Below we create an equilibrium problem with our
prescribed equilibrium conditions for *temperature*, *pressure*, and *amounts
of elements* that are consistent with our intention of calculating the
solubility of |CO2| at 60 °C and 100 bar in a 1 molal NaCl brine.

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step 5
    :end-before: Step 6

.. note::
    Did you pay attention we said prescribed equilibrium conditions for the
    *amounts of elements*? Since we actually provided the amounts of substances
    |H2O|, NaCl, and |CO2| in the above code, this statement seems a little bit
    confusing at least. Here is what happens behind the scenes: Reaktoro parses
    these chemical formulas and determines the elements and their coefficients.
    Once this is done, the amount of each element stored inside the object
    ``problem`` is incremented according to the given amount of substance and
    its coefficient in the formula. The amounts of elements you provide are
    then used as constraints for the Gibbs energy minimization calculation when
    computing the state of chemical equilibrium (i.e., when we try to find the
    amounts of all species in the system that correspond to a state of minimum
    Gibbs energy and at the same time satisfying the *element amounts
    constraints*).

.. danger::
    Now that you know that an equivalent chemical equilibrium problem could be defined with:

    .. code-block:: python

        problem.add('H' 111.0, 'mol')
        problem.add('O' 75.5, 'mol')
        problem.add('Na' 1.0, 'mol')
        problem.add('Cl' 1.0, 'mol')
        problem.add('C' 10.0, 'mol')

    assuming that 1 kg of |H2O| is roughly 55.5 mol, you might want to
    *adventure* in manually specifying different values for the amounts of
    elements. Just be extra careful with the values you provide as this could
    accidentally result in an *infeasible* chemical equilibrium state.

    Here is a simple example of element amount conditions that result in an
    infeasible equilibrium state. Consider a chemical system containing only a
    gaseous phase with gases |H2Og| and |CO2g|. Find non-negative amounts for
    |H2Og| and |CO2g| when the given amounts of elements are: 2 mol of H, 1 mol
    of O, and 1 mol of C.

.. note::
    Please note that we are not condemning the input form shown above in terms
    of element amounts, but only telling you to be attentive with the values
    you input. If you are using Reaktoro as a chemical reaction solver in a
    reactive transport simulator, for example, you'll most likely need to work
    directly with given amounts of elements, which shows that this input form
    is required in certain cases. For such time-dependent modeling problems, you
    often only need to ensure that the initial conditions for elements amounts
    result in feasible initial species amounts.

.. tip::
    The substance formulas given in the ``EquilibriumProblem.add`` method
    **can, but do not need to**, correspond to names of chemical species in the
    thermodynamic database. Even unusual, if not strange, substance formulas,
    such as HCl3(NaO)4C13, would be understood by that method. *We do not
    promise, however, that you will obtain a feasible chemical equilibrium
    state with unrealistic conditions!*

.. note::
    In Reaktoro, the word *element* is used as a synonym of components of
    chemical species, and not necessarily chemical elements. Electric charge
    is, for example, an element, even though it is not a not a chemical
    element. Thus, we say the ionic species CO\ :sub:`3`\ :sup:`2-` is composed
    of elements C, O, and Z, with coefficients 1, 3, and -2 respectively, where
    Z is the symbol we use to denote the electric charge element.

.. warning::
    Prefer the use of neutral substances when using the method
    ``EquilibriumProblem.add``, unless you definitely need to add a
    charged, ionic species in the recipe. The following code:

    .. code-block:: python

        problem.add('H+' 0.1, 'mmol')

    will increment not only the amount of element H by 0.1 mmol, but also the
    electric charge element Z. As a result, the composition of the aqueous
    phase at equilibrium will not be electrically neutral, which might not be
    your intention. The following:

    .. code-block:: python

        problem.add('H+' 0.1, 'mmol')
        problem.add('Cl-' 0.1, 'mmol')

    would result in the amount of element Z equal to zero.


Calculating the chemical equilibrium state
------------------------------------------

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step 6
    :end-before: Step 7

In this step, we use the ``equilibrate`` function to calculate the chemical
equilibrium state of the system with the given equilibrium conditions stored in
the object ``problem``. For this calculation, Reaktoro uses an efficient
**Gibbs energy minimization** computation to determine the species amounts that
correspond to a state of minimum Gibbs energy in the system, while satisfying
the prescribed amount conditions for temperature, pressure, and element
amounts. The result is stored in the object ``state``, of class
``ChemicalState``.

.. attention::

    In the vast majority of cases, you'll only have one object of class
    ``ChemicalSystem`` in your code and one or more objects of class
    ``ChemicalState`` describing different states of your chemical system!
    Reaktoro differentiates these two independent concepts: *chemical system
    definition* and *chemical system state*.

.. tip::

    The method ``equilibrate`` is a convenient function in Reaktoro. Consider
    using the class ``EquilibriumSolver`` for more advanced requirements. For
    example, if you have to perform many equilibrium calculations in sequence.
    The ``equilibrate`` method has a computational overhead because every call
    creates a new object of class ``EquilibriumSolver``. Preferably, this
    object should be created only once, and then used subsequently for all
    other equilibrium calculations. Here is a demonstration:

    .. code-block:: python

        solver = EquilibriumSolver(system) # Our chemical equilibrium solver
        state = ChemicalState(system)  # Our chemical state

        solver.solve(state, problem)  # Initial equilibrium calculation
        state.output('state-initial.txt')  # Output the initial equilibrium state

        problem.add('NaCl', 0.1, 'mol')  # Increment the amount of NaCl

        solver.solve(state, problem)  # Subsequent equilibrium calculation
        state.output('state-modified.txt')  # Output the modified equilibrium state




Outputting the calculated chemical state to a file
--------------------------------------------------

We have performed our chemical equilibrium calculation and now we want to
inspect the computed compositional state and its thermodynamic properties,
which can be done by outputting the chemical state to a file, here named
``result.txt``.

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step 7
    :end-before: Step 8

Printing the amounts of some aqueous species
--------------------------------------------

Here is just a small demonstration of getting species amount information from a
``ChemicalState`` object using the method ``ChemicalState.speciesAmount`` to
extract the amount of a few chemical species. Please inspect the API of
`ChemicalState`_ class to learn more about its methods.

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step 8
    :end-before: Step 9

Printing the amounts of element C in both aqueous and gaseous phases
--------------------------------------------------------------------

Finally, we print the amounts of element C in both aqueous and gaseous
phases:

.. literalinclude:: ../../../../demos/python/demo-equilibrium-co2-solubility-nacl-brine.py
    :start-at: Step 9

In this specific case in which there were no initial element C in the aqueous
phase, the value corresponding to the amount of element C in the aqueous phase
is our solubility of CO2 in the NaCl brine with the previously prescribed
conditions.

.. tip::
    If we had used 2 kg of |H2O|, we would have needed to divide the calculated
    amount of element C in the aqueous phase by 2 to obtain the solubility in
    molal (mol per kg of |H2O|).

.. _ChemicalState: https://reaktoro.org/cpp/classReaktoro_1_1ChemicalState.html