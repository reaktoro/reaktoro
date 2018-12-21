Tutorials: Chemical Kinetics
============================

In many cases, chemical equilibrium states are not sufficient in our analysis,
and chemical kinetics is then imperative for proper understanding of process.
Reaktoro can perform chemical kinetics simulations combined with chemical
equilibrium (i.e., part of the chemical system evolves under kinetics, while
the other is continuously in equilibrium at all times). This mode of
calculation is particularly useful for simulating chemically reaction systems
in which some reactions have rates that are many orders of magnitude higher
than others (and thus can be assumed in instantaneous equilibrium at any time).
Geochemical systems involving aqueous species and minerals are examples


Chemical equilibrium calculations are essential for many chemical reaction
modeling problems. Below is a list of tutorials demonstrating how Reaktoro's
chemical equilibrium modeling capabilities can be applied to solve different
problems.


.. toctree::
   :maxdepth: 1
   :glob:

   co2-solubility-nacl-brine
