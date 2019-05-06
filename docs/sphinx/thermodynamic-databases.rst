Thermodynamic Databases
=======================

Thermodynamic databases allow us to define and model a chemically reactive
system by providing the means for the computation of necessary thermodynamic
properties (e.g., standard Gibbs energies of species, equilibrium constants of
reactions). In such databases, we find a collection of chemical species and/or
reactions and their accompanying data, which include substance's name and
chemical formula, reaction equations, thermodynamic data and/or model
parameters.

.. attention::

    There are many thermodynamic databases available in the literature and they
    are in general very different from each other. At the moment, there is no
    standard on how thermodynamic databases should be formatted. Some databases
    are based on chemical reactions and their equilibrium constants (e.g.,
    PHREEQC databases), while others are based on substances and their model
    parameters for evaluation of their standard thermodynamic properties at
    temperature and pressure of interest (e.g., SUPCRT92 databases).

Reaktoro currently supports the following thermodynamic databases:

* SUPCRT92 databases;
* PHREEQC databases; and
* GEMS databases.

SUPCRT92 Databases
------------------

.. |supcrt98| replace:: :download:`supcrt98.xml <../../databases/supcrt/supcrt98.xml>`
.. |supcrt07| replace:: :download:`supcrt07.xml <../../databases/supcrt/supcrt07.xml>`
.. |supcrt98-organics| replace:: :download:`supcrt98-organics.xml <../../databases/supcrt/supcrt98-organics.xml>`
.. |supcrt07-organics| replace:: :download:`supcrt07-organics.xml <../../databases/supcrt/supcrt07-organics.xml>`
.. |slop98| replace:: :download:`slop98.dat <../../databases/supcrt/slop98.dat>`
.. |slop07| replace:: :download:`slop07.dat <../../databases/supcrt/slop07.dat>`

.. _supcrt98.xml: :download:`<../../databases/supcrt/supcrt98.xml>`
.. _slop98.dat: :download:`<../../databases/supcrt/slop98.dat>`
.. _slop07.dat: :download:`<../../databases/supcrt/slop07.dat>`
.. _supcrt07.xml: :download:`<../../databases/supcrt/supcrt07.xml>`
.. _supcrt98-organics.xml: databases/supcrt/supcrt98-organics.xml
.. _supcrt07-organics.xml: databases/supcrt/supcrt07-organics.xml

.. sidebar:: SUPCRT92 Database Files

    | |supcrt98|
    | |supcrt98-organics|
    | |supcrt07|
    | |supcrt07-organics|

The SUPCRT92 thermodynamic databases supported in Reaktoro are presented next.
They contain parameters for the calculation of standard thermodynamic
properties of aqueous species, gases, and minerals for temperatures 0-1000 °C
and pressures 1-5000 bar. The standard properties of aqueous species are
calculated using the revised *Helgeson-Kirkham-Flowers* (HKF) equations of
state and, for the gases and minerals, a thermodynamic model based on
Maier--Kelly heat capacity polynomial equation.

.. note::

    The thermodynamic databases ``supcrt98.xml`` and ``supcrt07.xml``, in XML
    format, were derived, respectively, from the original SUPCRT92 database
    files |slop98| and |slop07|. In the process, **all organic aqueous species
    were removed**! If you need them in your modeling problem, you should then
    use instead ``supcrt98-organics.xml`` and ``supcrt07-organics.xml``.

.. tip::

    If your problem requires an aqueous phase without organic species and you
    are using an automatic initialization scheme for its construction (e.g.,
    creating an aqueous phase with all species in the database whose elements
    are H, O, C, or Ca), then make sure you are using one of the SUPCRT92
    databases **without organic species**! Otherwise, you might end up with an
    aqueous phase containing an extremely long list of organic species that
    will only serve to decrease the performance of the calculations.

.. note::

    The equation of state of Wagner and Pruss (2002) is used to calculate the
    thermodynamic properties of water and its temperature and pressure
    derivatives.

PHREEQC Databases
-----------------

.. |Amm| replace:: :download:`Amm.dat <../../databases/phreeqc/Amm.dat>`
.. |frezchem| replace:: :download:`frezchem.dat <../../databases/phreeqc/frezchem.dat>`
.. |iso| replace:: :download:`iso.dat <../../databases/phreeqc/iso.dat>`
.. |llnl| replace:: :download:`llnl.dat <../../databases/phreeqc/llnl.dat>`
.. |minteq| replace:: :download:`minteq.dat <../../databases/phreeqc/minteq.dat>`
.. |minteq.v4| replace:: :download:`minteq.v4.dat <../../databases/phreeqc/minteq.v4.dat>`
.. |phreeqc| replace:: :download:`phreeqc.dat <../../databases/phreeqc/phreeqc.dat>`
.. |pitzer| replace:: :download:`pitzer.dat <../../databases/phreeqc/pitzer.dat>`
.. |sit| replace:: :download:`sit.dat <../../databases/phreeqc/sit.dat>`
.. |wateq4f| replace:: :download:`wateq4f.dat <../../databases/phreeqc/wateq4f.dat>`

Reaktoro can use PHREEQC_ as a :ref:`thermodynamic backend<Thermodynamic
Backends>`, which permits us to take advantage of the rich collection of
PHREEQC thermodynamic databases that are listed next.

.. sidebar:: PHREEQC Database Files

    | |Amm|
    | |frezchem|
    | |iso|
    | |llnl|
    | |minteq|
    | |minteq.v4|
    | |phreeqc|
    | |pitzer|
    | |sit|
    | |wateq4f|

.. _PHREEQC: https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/

GEMS Databases
--------------

Reaktoro can also use GEMS as a :ref:`thermodynamic backend<Thermodynamic
Backends>` and take advantage of its databases.

.. todo::

    Write about GEMS databases.
