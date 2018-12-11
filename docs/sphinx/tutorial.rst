Tutorial
========


.. attention::

    This tutorial is under construction. It should be finished within a few
    days.


Single-phase chemical equilibrium calculation
---------------------------------------------

We present below a Python script that performs a single-phase chemical
equilibrium calculation using Reaktoro. The single phase in the chemical system
is an *aqueous phase* with its chemical species automatically collected from a
given thermodynamic database. We are interested in calculating the amounts of
every aqueous species after we let 1 kg of H2O, 0.5 mol of NaCl, and 1 mol of
CO2 to react to a state of chemical equilibrium under constant temperature and
pressure of 60 °C and 100 bar respectively.

.. literalinclude:: examples/equilibrium-single-phase.py

We explain next every step in this chemical reaction modeling exercise.

Step 1
^^^^^^

.. literalinclude:: examples/equilibrium-single-phase.py
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

Step 2
^^^^^^

.. literalinclude:: examples/equilibrium-single-phase.py
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





Step 3
^^^^^^

.. literalinclude:: examples/equilibrium-single-phase.py
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

        editor.addAqueousPhase(['H2O(l)', 'H+', 'OH-', 'Na+', 'Cl-', 'HCO3-', 'CO3--', 'CO2(aq)'])

.. caution::

    If you manually specify the chemical species in a phase, you need to make
    sure that they exist in the thermodynamic database with **exact same
    name**! Replacing ``'CO2(aq)'`` above with ``'CO2'`` will cause an error if
    the database has no aqueous species with such name.

Step 4
^^^^^^

.. literalinclude:: examples/equilibrium-single-phase.py
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
of class ``ChemicalSystem`` is almost always needed to initialize those
classes! This is because without the information about the chemical system and
the means to evaluate thermodynamic properties, no chemical reaction
calculation, either equilibrium or kinetics, can be performed.

Step 5
^^^^^^

.. literalinclude:: examples/equilibrium-single-phase.py
    :start-at: Step 5
    :end-before: Step 6

We have now defined and constructed our chemical system of interest, enabling
us to move on to a next step in Reaktoro's modeling workflow: formulating and
solving interesting chemical reaction problems. In this tutorial we are
interested in computing the chemical equilibrium state of an aqueous phase for
given conditions of temperature and pressure. In addition to these, we also
provide the conditions for amounts of chemical elements at equilibrium. You can
observe the

In this step we use the class ``EquilibriumProblem`` to define a **chemical
equilibrium problem** in which **temperature**, **pressure**, and **element
amounts** are given and we seek the species amounts that correspond to a state
of chemical equilibrium.



(obtained from the , which consists of using the class
``EquilibriumProblem``.


We then specify:

- the **temperature** and **pressure**  at chemical equilibrium (60 °C and 100 bar); and
- the amount/mass of substances that should be reacted
together for which we are interested in the chemical equilibrium state (react
1.0 kg of H2O, 1.0 mol of CO2, and 0.5 mol of NaCl).

.. note::

    The chemical equilibrium state resulting from the reactions of given
    substances (in this case H2O, CO2, and NaCl) is not calculated here by
    following every compositional changes over time (chemical kinetic
    evolution) until all chemical species do not react anymore (i.e., the net
    rate of reactions are zero). Instead, Reaktoro performs an efficient
    **Gibbs energy minimization** calculation to determine the state of
    chemical equilibrium.

We note that the actual input conditions for chemical equilibrium calculations
are given in terms of *amounts of chemical elements*. Thus, once a substance and
its amount/mass are given, Reaktoro converts this input into equivalent amounts
of chemical elements. For example, 1 mol of CO2 is converted to 1 mol of C and 2
mol of O. That is why any substance formula can be used in the
``EquilibriumProblem.add`` method --- in the end, this formula is parsed, and
the chemical elements and their number of atoms determined (e.g., the string
``'CaCO3'`` is parsed into ``[('Ca', 1), ('C', 1), ('O', 3)]``). This is then
used to compute the amount of each element from the given amount of substance
(e.g., 1 mol of CaCO3 resulting in 1 mol of C, 1 mol of Ca, and 3
    mol of O).

.. tip::

    The substance formulas specified in the ``add`` method of class
    ``EquilibriumProblem`` **do not need** to correspond to names of chemical
    species in the thermodynamic database. Even unusual, if not strange,
    substance formulas, such as HCl3(NaO)4C13, would be understood by that
    method. *We do not promise, however, that you will obtain a feasible
    chemical equilibrium state with unrealistic conditions!*

Thus, the following:

.. code-block:: python

    problem.add('CO2', 1.0, 'mol')

is equivalent to:

.. code-block:: python

    problem.add('C', 1.0, 'mol')
    problem.add('O', 2.0, 'mol')


In the above tip, we learn that Reaktoro converts amounts/masses of substances
into amounts of chemical elements. This is done observing the stoichiometry of
the elements in each substance. The following table shows the amounts of
elements resulting from the combination of 1 kg of H2O, 1 mol of CO2, and 0.5
mol of NaCl, assuming that 55 mol of H2O is roughly 1 kg of H2O.

.. table:: Amounts of chemical elements (approximately) obtained from the
    recipe 1.0 kg of H2O, 1.0 mol of CO2, and 0.5 mol of NaCl.
    :widths: 50 50
    :align: center

    ======= ============
    Element Amount (mol)
    ======= ============
    H       110.0
    O       57.0
    Na      0.5
    Cl      0.5
    C       1.0
    ======= ============

.. danger::

    Now that you know that an equivalent chemical equilibrium problem could be defined with:

    .. code-block:: python

        problem.add('H' 110.0, 'mol')
        problem.add('O' 57.0, 'mol')
        problem.add('Na' 0.5, 'mol')
        problem.add('Cl' 0.5, 'mol')
        problem.add('C' 1.0, 'mol')

    you might want to "adventure" in manually specifying different values for the amounts of elements.

    Note, however, that this input form has to be done carefully to ensure that the given element amounts result in a feasible chemical equilibrium state. Here is an example of conditions that result in unfeasible equilibrium states for a specific chemical system: calculate equilibrium of gases H2O and CO2 where the amount of H is 2 mol, O is 1 mol, and C is 1 mol.

Step 6
^^^^^^

.. literalinclude:: examples/equilibrium-single-phase.py
    :start-at: Step 6
    :end-before: Step 7

Step 7
^^^^^^

.. literalinclude:: examples/equilibrium-single-phase.py
    :start-at: Step 7
    :end-before: Step 8

Step 8
^^^^^^

.. literalinclude:: examples/equilibrium-single-phase.py
    :start-at: Step 8
