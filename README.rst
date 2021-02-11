========
Reaktoro
========

.. image:: https://github.com/reaktoro/reaktoro/workflows/linux/badge.svg?branch=master
    :alt: Linux Build
    :scale: 100%
    :target: https://github.com/reaktoro/reaktoro/actions?query=workflow%3Alinux

.. image:: https://github.com/reaktoro/reaktoro/workflows/osx/badge.svg?branch=master
    :alt: OSX Build
    :scale: 100%
    :target: https://github.com/reaktoro/reaktoro/actions?query=workflow%3Aosx

.. image:: https://github.com/reaktoro/reaktoro/workflows/windows/badge.svg?branch=master
    :alt: Windows Build
    :scale: 100%
    :target: https://github.com/reaktoro/reaktoro/actions?query=workflow%3Awindows


Reaktoro is a unified framework for modeling chemically reactive systems. It provides methods for chemical equilibrium and kinetic calculations for multiphase systems. Reaktoro is mainly developed in C++ for performance reasons. A Python interface is available for a more convenient and simpler use. Currently, Reaktoro can interface with two widely used geochemical software: `PHREEQC <http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/>`_ and `GEMS <http://gems.web.psi.ch/>`_.

Quick example
=============

Here is a simple C++ code using Reaktoro to perform a multiphase chemical
equilibrium calculation:

.. code-block:: cpp

    #include <Reaktoro/Reaktoro.hpp>

    using namespace Reaktoro;

    int main()
    {
        ChemicalEditor editor;
        editor.addAqueousPhase("H2O NaCl CaCO3 CO2");
        editor.addGaseousPhase("CO2(g)");
        editor.addMineralPhase("Calcite");

        ChemicalSystem system(editor);

        EquilibriumProblem problem(system);
        problem.add("H2O", 1, "kg");
        problem.add("CO2", 1, "mol");
        problem.add("NaCl", 0.7, "mol");
        problem.add("CaCO3", 1, "g");

        ChemicalState state = equilibrate(problem);

        state.output("result.txt");
    }


This calculation could also be performed using Reaktoro's Python interface:

.. code-block:: python

    from reaktoro import *

    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O NaCl CaCO3 CO2")
    editor.addGaseousPhase("CO2(g)")
    editor.addMineralPhase("Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("CO2", 1, "mol")
    problem.add("NaCl", 0.7, "mol")
    problem.add("CaCO3", 1, "g")

    state = equilibrate(problem)

    state.output("result.txt")


Installation and Tutorials
==========================

For installation instructions, tutorials, and list of publications related to
this project, please access `reaktoro.org <http://www.reaktoro.org>`_. This web
site describes how to download and install Reaktoro, and demonstrate some basic
usage.

FAQ
===

How do I ask a question about Reaktoro?
---------------------------------------

If you have questions about using or installing Reaktoro, please go to
`Reaktoro's GitHub Issues`_ and let us know. Please select the **question**
label on the right side of the issue pages. We'll do our best to answer your
question as soon as possible.


How can I report a bug?
-----------------------

You got a bug and this is frustrating, we understand you. But don't worry —
we'll be happy to fix it for you (*provided it is indeed a bug!*).

Before you report a bug, please check first if someone else has already
reported the same issue. If not, go to `Reaktoro's GitHub Issues`_ and enter a
*descriptive title* and *write your issue with enough details*. Please select
the label **bug** on the right side of the page.

Please provide a `Minimum Reproducible Example`_? Please provide such an
example so that we can be more efficient in identifying the bug and fixing it
for you.

Have you heard about `Markdown`_? Please use Markdown syntax when reporting
your issues.

How can I contribute to Reaktoro?
---------------------------------

First, thanks for your interest in contributing to Reaktoro! You can do so in
many ways, from reporting bugs and writing tutorials to helping us with code
development. You might also consider **financially supporting Reaktoro's
development** by helping us extending the development team if you plan to make
Reaktoro an essential software component in your company or academic group.

Read more on how to contribute to Reaktoro `here <CONTRIBUTING.rst>`__.

Contributors
============

You can see the list of awesome people who has contributed code to Reaktoro in
the `contributors page
<https://github.com/reaktoro/Reaktoro/graphs/contributors>`__.

We would love to have you as a contributor too, see `CONTRIBUTING
<CONTRIBUTING.rst>`__ for more information.

Developing Quick-Start
======================

In order to start developing, you'll need to build Reaktoro from sources. There
are two ways: install the dependencies manually, as described `here
<http://www.reaktoro.org/installation.html>`_, or using Conda. `Conda
<https://conda.io/docs/>`_ is a tool for managing packages, dependencies and
environments for multiple languages, including Python and C++, and supporting
multiple platforms: Windows, Linux and macOS. In order to start developing
Reaktoro using Conda, these are the steps:

#. Install Miniconda, pick the 64-bit installer that uses the latest Python version from: `conda.io/miniconda.html <https://conda.io/miniconda.html>`_.
#. Add ``conda-forge`` as a channel: ``conda config --append channels conda-forge``
#. Install ``conda-devenv``: ``conda install -n base conda-devenv``
#. Create an environment for Reaktoro, from the repository root directory: ``conda devenv``
#. Activate the environment: ``source activate reaktoro`` from Linux/macOS or ``activate reaktoro`` from Windows
#. Create a ``build`` directory and call ``cmake`` from it (for now check the `.travis.yml` file for an example on CMake parameters), OR, on Windows, call the ``inv msvc`` task to generate a project under ``build\msvc`` directory, open it in the IDE and build the ``INSTALL`` project. (``inv`` is short for ``invoke``, from the `Invoke <https://www.pyinvoke.org/>`_ tool.)


License
=======

LGPL v2.1

Copyright (C) 2014-2018 Allan Leal

Reaktoro is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option) any
later version.

Reaktoro is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.


.. _Reaktoro's GitHub Issues: https://github.com/reaktoro/Reaktoro/issues/new
.. _Minimum Reproducible Example: https://stackoverflow.com/help/mcve>
.. _Markdown: https://guides.github.com/features/mastering-markdown/

__ `Reaktoro's GitHub Issues`_
