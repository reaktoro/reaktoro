========
Reaktoro
========

.. image:: https://travis-ci.org/reaktoro/Reaktoro.svg?branch=master
    :target: https://travis-ci.org/reaktoro/Reaktoro

.. image:: https://ci.appveyor.com/api/projects/status/github/reaktoro/Reaktoro?branch=master&svg=true
    :target: https://ci.appveyor.com/project/reaktoro/Reaktoro


Reaktoro is a unified framework for modeling chemically reactive systems. It provides methods for chemical equilibrium and kinetic calculations for multiphase systems. Reaktoro is mainly developed in C++ for performance reasons. A Python interface is available for a more convenient and simpler use. Currently, Reaktoro can interface with two widely used geochemical software: `PHREEQC <http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/>`_ and `GEMS <http://gems.web.psi.ch/>`_.

Quick example
=============

Here is a simple C++ code using Reaktoro to perform a multiphase chemical equilibrium calculation:

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

For installation instructions, tutorials, and list of publications related to this project, please access `reaktoro.org <http://www.reaktoro.org>`_. This web site describes how to download and install Reaktoro, and demonstrate some basic usage.

FAQ
===
Do you have a question about Reaktoro?
--------------------------------------
Please open an issue at [Reaktoro's GitHub Issues][github-issues] using the **question** label on the right side of the issue pages. That one of ours developers will answer as soon as possible.

Did you find a bug?
-------------------------
You got a bug and this is frustrating, we understand you. But don't worry - we'll be happy to fix it for you (*provided it is indeed a bug!*). 
Please open an issue at the [Reaktoro's GitHub Issues][github-issues] using the **bug** label ([read more]()).

Want to contribuite with Reaktoro?
----------------------------------
First, thanks for your interest in contributing with Reaktoro! We love feedback in all forms: issues, comments, PRs, etc.
[Click here](CONTRIBUTING.md) to find out more about our workflow.

License
=======

LGPL v2.1

Copyright (C) 2014-2018 Allan Leal

Reaktoro is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

Reaktoro is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
