Tutorial
========

Single-phase chemical equilibrium calculation
---------------------------------------------

We present below a python script that performs a single-phase chemical
equilibrium calculation using Reaktoro. The single phase in the chemical system
is an *aqueous phase* with chemical species automatically collected from a
given thermodynamic database We then explain every step in this modeling
exercise.

.. literalinclude:: examples/equilibrium-single-phase.py

Step 1
^^^^^^

Using Reaktoro in python requires that the ``reaktoro`` python
package be imported.

.. note::

    To simplify the tutorials, ``from reaktoro import *`` is used, which
    imports all components of the **reaktoro** package. For your applications,
    you might consider using ``import reaktoro as rkt`` instead, and then refer
    to Reaktoro's methods and classes with ``rkt::ChemicalSystem``,
    ``rkt::equilibrate``, and so forth.

Step 2
^^^^^^

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

Step 3
^^^^^^

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

Step 4
^^^^^^

Once we have our chemical system defined and constructed, we can now move on to
defining interesting chemical reaction problems to solve. In this step, we
define a chemical equilibrium problem using the class ``EquilibriumProblem``.
We then specify *temperature* and *pressure conditions* at equilibrium (60 °C
and 100 bar), as well as the amount/mass of substances that should be reacted
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

.. tip::

    The chemical formulas of substances specified in the ``add`` method of
    class ``EquilibriumProblem`` **do not need** to correspond to names of
    chemical species in the thermodynamic database. Even non-existent substance
    formulas, such as HCl3(NaO)4C13, would be understood by the
    ``EquilibriumProblem.add`` method. *We do not promise, however, that you will obtain a feasible chemical
    equilibrium state with unrealistic conditions!*

    The reason for this flexibility is because only the *chemical elements* in
    these substance formulas (H2O, CO2, NaCl), and their *number of atoms*, are
    actually required to specify the mass conservation conditions of a chemical
    equilibrium problem.     , which are automatically obtained from Reaktoro's elemental
    formula parser. For example, the string ``'H2O'`` when parsed produces a
    list ``[('H', 2), ('O', 1))]``. Assuming 55 mol of H2O roughly corresponds
    to 1 kg of this substance, this input condition will be transformed into
    **110 mol of H** and **55 mol of O**.



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


Step 5
^^^^^^

Step 6
^^^^^^

Step 7
^^^^^^

