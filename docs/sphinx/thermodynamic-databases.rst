Thermodynamic Databases
=======================

Reaktoro currently supports the following thermodynamic databases:

* SUPCRT92 databases; and
* PHREEQC databases.

SUPCRT92 Databases
------------------

.. |supcrt98| replace:: :download:`supcrt98.xml <../../databases/supcrt/supcrt98.xml>`
.. |supcrt07| replace:: :download:`supcrt07.xml <../../databases/supcrt/supcrt07.xml>`
.. |supcrt98-organics| replace:: :download:`supcrt98-organics.xml <../../databases/supcrt/supcrt98-organics.xml>`
.. |supcrt07-organics| replace:: :download:`supcrt07-organics.xml <../../databases/supcrt/supcrt07-organics.xml>`
.. |slop98| replace:: :download:`slop98.dat <../../databases/supcrt/slop98.dat>`
.. |slop07| replace:: :download:`slop07.dat <../../databases/supcrt/slop07.dat>`

.. sidebar:: SUPCRT92 Database Files

    | |supcrt98|
    | |supcrt07|
    | |supcrt98-organics|
    | |supcrt07-organics|

SUPCRT92 databases contain parameters for the calculation of standard
thermodynamic properties of aqueous species, gases, and minerals for
temperatures 0-1000 °C and pressures 1-5000 bar. The standard properties of
aqueous species are calculated using the revised *Helgeson-Kirkham-Flowers*
(HKF) equations of state and, for the gases and minerals, a thermodynamic model
based on Maier--Kelly heat capacity polynomial equation.

.. csv-table:: A brief description of SUPCRT92 databases in Reaktoro.
    :header: "Reaktoro Database", "Brief Description"

    |supcrt98|,          "Derived from the SUPCRT92 database |slop98|, **without** organic species."
    |supcrt07|,          "Derived from the SUPCRT92 database |slop07|, **without** organic species."
    |supcrt98-organics|, "Derived from the SUPCRT92 database |slop98|, **with** organic species."
    |supcrt07-organics|, "Derived from the SUPCRT92 database |slop07|, **with** organic species."


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

Reaktoro can use PHREEQC as a thermodynamic back-end. This feature allows
Reaktoro users to use not only the thermodynamic databases of PHREEQC, but also
its activity models. Thus, by using PHREEQC as a thermodynamic back-end in
Reaktoro, users are able to take advantage of its rich collection of databases,
with the standard ones listed next.

.. _supcrt98.xml: :download:`<../../databases/supcrt/supcrt98.xml>`
.. _slop98.dat: :download:`<../../databases/supcrt/slop98.dat>`
.. _slop07.dat: :download:`<../../databases/supcrt/slop07.dat>`
.. _supcrt07.xml: :download:`<../../databases/supcrt/supcrt07.xml>`
.. _supcrt98-organics.xml: databases/supcrt/supcrt98-organics.xml
.. _supcrt07-organics.xml: databases/supcrt/supcrt07-organics.xml
