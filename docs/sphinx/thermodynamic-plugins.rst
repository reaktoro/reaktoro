Thermodynamic Plugins
=====================

There are excellent open-source software for modeling chemical systems. For
geochemical modeling in particular, two widely used software are PHREEQC_ and
GEMS_.

.. figure:: img/phreeqc-logo.svg
    :width: 50%
    :align: center
    :target: PHREEQC_

    PHREEQC is a computer program for speciation, batch-reaction,
    one-dimensional transport, and inverse geochemical calculations developed
    by USGS, United States.

.. figure:: img/gems-logo.png
    :width: 50%
    :align: center
    :target: GEMS_

    GEMS is a Gibbs energy minimization software for geochemical modeling
    developed at Paul Scherrer Institute, Switzerland.

Because many scientists and academics rely on specific thermodynamic models
(e.g., models for activity coefficients, fugacity coefficients) and databases
(e.g., PSI/NAGRA, SUPCRT92, PHREEQC databases) supported by these two software,
Reaktoro has been developed to permit the use of PHREEQC and GEMS as
thermodynamic plugins. This allows you to take advantage of the best there is
in these classical modeling codes together with fast and robust numerical
algorithms implemented in Reaktoro for intense chemical reaction calculations.

.. note::

    Neither PHREEQC nor GEMS are actually used for chemical reaction
    calculations (i.e., for chemical equilibrium/kinetics calculation) in
    Reaktoro. These plugins act as providers of thermodynamic properties of
    species, phases, and reactions (e.g., activity coefficients, activities,
    standard chemical potentials, phase molar volume and enthalpy, equilibrium
    constant of reaction). Whenever Reaktoro's fast algorithms need those
    properties, these plugins are invoked to retrieve them in a very efficient
    way by directly using their API.

PHREEQC Plugin
--------------

The code below demonstrates the combined use of Reaktoro and PHREEQC to perform
a chemical equilibrium calculation in which PHREEQC thermodynamic data and
activity models are used together with Reaktoro's Gibbs energy minimization
algorithm.

.. literalinclude:: ../../demos/python/demo-phreeqc-ex1.py
    :start-at: from
    :language: python

To execute the above code, assuming you have already installed Reaktoro, you
can download the Python script
:download:`demo-phreeqc-ex1.py<../../demos/python/demo-phreeqc-ex1.py>` and the
database file :download:`phreeqc.dat<../../databases/phreeqc/phreeqc.dat>` to
the same directory and execute in the terminal:

.. code::

    python demo-phreeqc-ex1.py

.. code-block:: c++

    #include <Reaktoro/Reaktoro.hpp>
    using namespace Reaktoro;

    int main()
    {
        // Initialize a Phreeqc object with a PHREEQC database file such as
        // phreeqc.dat, pitzer.dat, llnl.dat, sit.dat, minteq.dat, etc.
        Phreeqc phreeqc;("some-phreeqc-database-file.dat");

        // Use the just created Phreeqc object to execute a PHREEQC script,
        phreeqc.execute("some-phreeqc-input-script.dat");

        // and then use it to construct a ChemicalSystem object containing all
        // phases and species selected by PHREEQC.
        ChemicalSystem system = phreeqc;

        // Create a ChemicalState object for the PHREEQC calculated chemical state,
        // which was produced when the Phreeqc object executed the script above.
        ChemicalState state = phreeqc.state(system);

        // Define an equilibrium problem that uses the calculated PHREEQC state as
        // a basis for the calculation of a new equilibrium state. The problem
        // below adds 0.1 moles of CO2 to the calculated PHREEQC state while
        // maintaining the same temperature and pressure specified in the script.
        EquilibriumProblem problem(system);
        problem.setTemperature(state.temperature());
        problem.setPressure(state.pressure());
        problem.add(state);
        problem.add("CO2", 0.1, "moles");

        // Calculate the new equilibrium state using the above equilibrium problem.
        equilibrate(state, problem);

        // Output the updated equilibrium state to a file.
        state.output("state.txt");
    }

GEMS Plugin
-----------

The code below briefly demonstrates how GEMS thermodynamic models and databases
can be used with Reaktoro's numerical methods.

.. code-block:: c++

    #include <Reaktoro/Reaktoro.hpp>
    using namespace Reaktoro;

    int main()
    {
        // Use an exported project file from GEMS to initialize a Gems object,
        Gems gems;("project.lst");

        // and then use it to construct the ChemicalSystem object.
        ChemicalSystem system = gems;

        // Create a ChemicalState object that contains the temperature, pressure,
        // and amounts of species stored in the exported GEMS file.
        ChemicalState state = gems.state(system);

        // Change the temperature of the chemical state,
        state.setTemperature(50, "celsius");

        // and then equilibrate the modified chemical state using Reaktoro's methods.
        equilibrate(state);

        // Output the updated equilibrium state to a file.
        state.output("state.txt");
    }

What about more thermodynamic plugins?
--------------------------------------

Are there other chemical reaction modeling software that you think could be
integrated with Reaktoro? Let us know by creating a new issue at `Reaktoro's
GitHub Issues`_.

.. attention::

    It would be great if you could contribute to expanding the list of supported
    Reaktoro's thermodynamic plugins. Contributions can be made in several
    forms, ranging from direct code contribution to financing a project in
    which one or more experts will accomplish this.


.. _PHREEQC: http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
.. _GEMS: http://gems.web.psi.ch/
.. _Reaktoro's GitHub Issues: https://github.com/reaktoro/reaktoro/issues/new
