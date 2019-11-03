Tutorials: Chemical Kinetics
============================

Chemical equilibrium calculations alone are sometimes not sufficient to
understand a chemically reactive process. This happens when we need to
understand how the composition of the chemical system changes with time as a
result of chemical reactions. For this, *chemical kinetics* is imperative.

Reaktoro can perform chemical kinetics calculations combined with chemical
equilibrium (i.e., part of the chemical system evolves under kinetics, while
the other is continuously in equilibrium at all times). This mode of
calculation is particularly useful for simulating chemically reactive systems
in which some reactions have rates that are many orders of magnitude higher
than others (and thus can be assumed in instantaneous equilibrium at any time).

Geochemical systems involving aqueous species and minerals are examples in
which this combined chemical kinetics and equilibrium approach is valuable.

.. toctree::
   :maxdepth: 1
   :glob:

   kineticpath-carbonates-co2
   kineticpath-calcite-hcl
