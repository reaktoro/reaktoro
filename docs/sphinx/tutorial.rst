Tutorial
========

Single-phase chemical equilibrium calculation
---------------------------------------------

We present below a python script that performs a single-phase chemical equilibrium calculation using Reaktoro and we subsequently explain every step in the modeling exercise.

.. code-block:: python

    # Step 1: Import the reaktoro python package
    from reaktoro import *

    # Step 2: Define your chemical system
    editor = ChemicalEditor()
    editor.addAqueousPhase("H O Na Cl C")

    # Step 3: Construct your chemical system
    system = ChemicalSystem(editor)

    # Step 4: Define the chemical equilibrium problem
    problem = EquilibriumProblem(system)
    problem.setTemperature(60, 'celsius')
    problem.setPressure(100, 'bar')
    problem.add('H2O', 1.0, 'kg')
    problem.add('NaCl', 0.5, 'mol')
    problem.add('CO2', 1.0, 'mol')

    # Step 5: Calculate the chemical equilibrium state
    state = equilibrate(problem)

    # Step 6: Output the calculated chemical state to a file
    state.output('result.txt')

    # Step 7: Print the mole amounts of some aqueous species
    print('Amount of CO2(aq): ', state.speciesAmount('CO2(aq)'))
    print('Amount of HCO3-: ', state.speciesAmount('HCO3-'))
    print('Amount of CO3--: ', state.speciesAmount('CO3--'))

Step 1
^^^^^^

Step 2
^^^^^^

Step 3
^^^^^^

Step 4
^^^^^^

Step 5
^^^^^^

Step 6
^^^^^^

Step 7
^^^^^^

