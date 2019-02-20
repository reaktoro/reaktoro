Thermodynamic Backends
======================

There are excellent and mature open-source software for modeling chemical
systems for which development and maintenance has been ongoing for several
years to a few decades. Wouldn't it be great if Reaktoro could use the best
these software can offer in terms of thermodynamic databases and activity
models as *thermodynamic backends*, while relying on Reaktoro's numerical
algorithms for the intense chemical reaction calculations? Well, it turns out
this is already supported for two widely used geochemical modeling codes:
PHREEQC_ and GEMS_.

.. figure:: img/phreeqc-logo.svg
    :figwidth: 80%
    :width: 50%
    :align: center
    :target: PHREEQC_

    PHREEQC is a computer program for speciation, batch-reaction,
    one-dimensional transport, and inverse geochemical calculations developed
    by USGS, United States.

.. figure:: img/gems-logo.png
    :figwidth: 80%
    :width: 50%
    :align: center
    :target: GEMS_

    GEMS is a Gibbs energy minimization software for geochemical modeling
    developed at Paul Scherrer Institute, Switzerland.

.. note::

    Neither PHREEQC nor GEMS are used in Reaktoro for solving the underlying
    mathematical problems for chemical equilibrium and kinetics. These backends
    act as providers of thermodynamic properties of species, phases, and
    reactions (e.g., activity coefficients, activities, standard chemical
    potentials, phase molar volume and enthalpy, equilibrium constant of
    reaction). Whenever Reaktoro's numerical algorithms need these properties,
    the thermodynamic backend is invoked to retrieve them in a very efficient
    way by directly using its API. It seems complicated, but rest assure, the
    usage is rather simple as you'll see below!

PHREEQC Backend
---------------

The code below demonstrates the combined use of Reaktoro and PHREEQC to perform
a chemical equilibrium calculation in which PHREEQC thermodynamic data and
activity models are used together with Reaktoro's Gibbs energy minimization
algorithm.

.. literalinclude:: ../../demos/python/demo-backends-phreeqc.py
    :start-at: from
    :language: python

Python and C++ files for this demo:
    | :download:`demo-backends-phreeqc.py<../../demos/python/demo-backends-phreeqc.py>`
    | :download:`demo-backends-phreeqc.cpp<../../demos/cpp/demo-backends-phreeqc.cpp>`

GEMS Backend
------------

Similarly, the code below briefly demonstrates how Reaktoro and GEMS can be
used together. You'll need first to prepare your chemical system definition
using GEM-Selektor, the graphical user interface of GEMS. In this step, you'll
be able to select which GEMS' supported thermodynamic database you want to use
as well as the activity models for each phase (aqueous, gaseous, solid
solutions). Next, export the GEMS project files to disk, and use it in Reaktoro
as shown below.

.. literalinclude:: ../../demos/python/demo-backends-gems.py
    :start-at: from
    :language: python

Python and C++ files for this demo:
    | :download:`demo-backends-gems.py<../../demos/python/demo-backends-gems.py>`
    | :download:`demo-backends-gems.cpp<../../demos/cpp/demo-backends-gems.cpp>`

What about more thermodynamic backends?
---------------------------------------

Are there other chemical reaction modeling software that you think could be
integrated with Reaktoro as *thermodynamic backends*? Let us know by creating a
new issue at `Reaktoro's GitHub Issues`_.

.. attention::
    It would be great if you could contribute to expanding the list of
    supported Reaktoro's thermodynamic backends. Contributions can be made in
    several forms, ranging from direct code contribution to financing a project
    in which one or more experts will implement this.

.. _PHREEQC: http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/
.. _GEMS: http://gems.web.psi.ch/
.. _Reaktoro's GitHub Issues: https://github.com/reaktoro/reaktoro/issues/new
