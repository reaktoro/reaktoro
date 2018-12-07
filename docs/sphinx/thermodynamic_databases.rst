Thermodynamic Databases
=======================

The currently supported thermodynamic databases in Reaktoro are listed
here.

SUPCRT92 Databases
------------------

current The table below shows the current available built-in databases
in Reaktoro and their description.

+-----------------------------------+-----------------------------------+
| Built-in Database                 | Description                       |
+===================================+===================================+
| `supcrt98.xml`_                   | Derived from the                  |
|                                   | SUPCRT92@sup{@cite Johnson1992}   |
|                                   | database file `slop98.dat`_,      |
|                                   | **without** organic species.      |
+-----------------------------------+-----------------------------------+
| `supcrt07.xml`_                   | Derived from the                  |
|                                   | SUPCRT92@sup{@cite Johnson1992}   |
|                                   | database file `slop07.dat`_,      |
|                                   | **without** organic species.      |
+-----------------------------------+-----------------------------------+
| `supcrt98-organics.xml`_          | Derived from the                  |
|                                   | SUPCRT92@sup{@cite Johnson1992}   |
|                                   | database file `slop98.dat`_,      |
|                                   | **with** organic species.         |
+-----------------------------------+-----------------------------------+
| `supcrt07-organics.xml`_          | Derived from the                  |
|                                   | SUPCRT92@sup{@cite Johnson1992}   |
|                                   | database file `slop07.dat`_,      |
|                                   | **with** organic species.         |
+-----------------------------------+-----------------------------------+

This file was produced by converting database has been produced for
Reaktoro is produced from the is the SUPCRT92@sup{@cite Johnson1992}
database which contains parameters for the revised
*Helgeson-Kirkham-Flowers* (HKF) equations of state for the calculation
of standard thermodynamic properties of hundreds of aqueous species as
well as Maier–Kelly coefficients for the calculation of standard
thermodynamic properties of gases and minerals at temperatures 0 to 1000
°C and pressures 1 to 5000 bar.

that should be found **at the same directory from where the application
is executed!** For convenience, Reaktoro also maintains a few built-in
databases, including ``supcrt98.xml``, whose files need not be found
along with the application. Thus, if no file ``supcrt98.xml`` is found,
the `Database`_ object will be initialized with the built-in database
instead.

PHREEQC Databases
-----------------

Reaktoro supports PHREEQC@sup{@cite Parkhurst2013} thermodynamic
databases.


.. _supcrt98.xml: databases/supcrt/supcrt98.xml
.. _slop98.dat: databases/supcrt/slop98.dat
.. _supcrt07.xml: databases/supcrt/supcrt07.xml
.. _slop07.dat: databases/supcrt/slop07.dat
.. _supcrt98-organics.xml: databases/supcrt/supcrt98-organics.xml
.. _supcrt07-organics.xml: databases/supcrt/supcrt07-organics.xml
.. _Database: @ref%20Reaktoro::Database