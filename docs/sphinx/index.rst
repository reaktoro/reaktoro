.. Reaktoro documentation master file, created by
   sphinx-quickstart on Tue Dec  4 08:56:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. Reaktoro -- a unified framework for modeling chemically reactive systems
.. ========================================================================

Reaktoro
========

.. figure:: img/reaktoro-header.svg
   :align: center
   :width: 100%

:emphasis:`Reaktoro is a computational framework developed in C++ and Python
that implements numerical methods for modeling chemically reactive processes
governed by either chemical equilibrium, chemical kinetics, or a combination of
both.`

.. note::
    This guide is under active development. You may find some pages incomplete
    at the moment. *17.12.2018*

.. only: not latex

.. toctree::
    :caption: Getting Started
    :maxdepth: 1

    overview
    installation

.. toctree::
    :caption: Tutorials
    :maxdepth: 1

    Chemical Equilibrium <tutorials/equilibrium/index>
    Chemical Kinetics <tutorials/kinetics/index>
    Chemical Transport <tutorials/transport/index>

.. toctree::
    :caption: General
    :maxdepth: 1

    thermodynamic-databases
    thermodynamic-backends

.. toctree::
    :caption: About
    :maxdepth: 1

    acknowledgements
    license
