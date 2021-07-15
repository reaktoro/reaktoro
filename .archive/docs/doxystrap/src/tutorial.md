# Tutorial {#tutorial}

In this tutorial, we'll cover the main classes and methods in %Reaktoro for performing multiphase chemical equilibrium and kinetics calculations.

## Introduction

Before you start reading the tutorials, it might be worth reading this introductory section containing essential information for the understanding of the underlying theory for modeling chemically reactive systems using either chemical equilibrium or chemical kinetics models. 

This introduction is not meant to be comprehensive. We recommend you to read the following article for a more in-depth discussion about key concepts in chemical equilibrium and chemical kinetics:

<div class="list-group">
    <a href="http://doi.org/10.1515/pac-2016-1107" class="list-group-item">
      <p class="list-group-item-text">Leal, A.M.M., Kulik, D.A., Smith, W.R., Saar, M.O. <b>(2017)</b>. *An overview of computational methods for chemical equilibrium and kinetic calculations for geochemical and reactive transport modeling*. **Pure and Applied Chemistry**, 89(5), 597–643.</p>
    </a>
</div>

### Chemical System, Phases, Species, and Elements

When modeling chemically reactive processes, one has to define a *chemical system*. The chemical system definition comprehends all possible phases that could exist for the modeling problem of interest. Examples of phases include aqueous, gaseous, liquid, and solid solutions, pure minerals, plasma. 

Which phases should one consider when defining the chemical system? Selecting an appropriate set of phases not always is a straightforward task. For example, when modeling the chemical reactions within an aqueous solution, either using chemical equilibrium or chemical kinetics models, it is clear that the chemical system should consider an *aqueous phase*. However, suppose this computational modeling is performed at various temperature conditions so that the following possible scenarios can exist:

* What if a decrease in temperature causes the aqueous solution to precipitate one or more mineral phases? 
* What if an increase in temperature causes gases to exsolve? 
* What if the temperature increase is so high that the aqueous solution fully evaporates to become a gaseous solution and, as a result, precipitates several minerals that were previously dissolved? 

The chemical calculations can only predict these phase appearance and disappearance behavior as long as all potentially existing phases are considered in the definition of your chemical system!

Each phase in a chemical system is composed by one or more *chemical species*. Examples of chemical species include:

* *ionic aqueous species*, e.g., Na@sup{+}(aq), Cl@sup{-}(aq), HCO@sub{3}@sup{-}(aq);
* *neutral aqueous species*, e.g., CO@sub{2}(aq), H@sub{2}O(l);
* *gases*, e.g., CO@sub{2}(g), H@sub{2}S(g), CH@sub{4}(g), N@sub{2}(g); 
* *minerals*, e.g., CaCO@sub{3}(s), SiO@sub{2}(s), Al@sub{2}Si@sub{2}O@sub{5}(OH)@sub{4}(s));
* *complex organic molecules*, e.g., proteins, lipids, carbohydrates, vitamins.

Each chemical species is composed by one or more *components*. Components can be *chemical elements* (e.g., H, O, C, Na, Cl, Ca, Si), *electrical charge* (Z), as well as a combination of chemical elements and electrical charge, commonly known as *primary species* (e.g., H@sup{+}(aq), H@sub{2}O(l), CO@sub{2}(aq), HCO@sub{3}@sup{-}(aq), Fe@sup{2+}(aq)).

@note The word component has many meanings in the scientific literature. It is sometimes used to actually denote the chemical species, rather than the entities that constitute them. To avoid any misunderstanding between species and components, we prefer the use of word *element*. Thus, an element can denote chemical elements, electrical charge, or even primary species. The chemical species are then said to be composed by one or more elements.

### Chemical Reactions: How to model them?

<!-- Chemical reactions can be *homogeneous reactions*, when the reacting species belong to the same phase, or *heterogeneous reactions*, when the reacting species belong to different phases. -->

Chemical reactions can be modeled by using either *chemical kinetics*, *chemical equilibrium*, or a combination of both. In all cases, an initial condition for the chemical state of the system (e.g., temperature, pressure, amounts of the chemical species) are required to calculate its final state of interest. When using chemical equilibrium for modeling the reactions, this final state corresponds to a state in which the system is in chemical equilibrium, with the species concentrations no longer experiencing any changes with time. When using chemical kinetics, or a combination of chemical kinetics and equilibrium, the final state of interest might or not be the one corresponding to chemical equilibrium.

A chemical kinetics model for the reactions is capable of tracing the amounts of the chemical species over time as they undergo a series of chemical reactions. Under certain conditions (e.g., in a system closed to mass inflow/outflow) these chemical reactions might tend to a state of equilibrium, i.e., the forward and reverse rates of the chemical reactions balance and the species amounts no longer change with time. It is possible, however, that the chemical system tend first to a *metastable equilibrium* state, as a result of energy barriers that prevent the reactions to continue towards the most stable equilibrium state. Furthermore, it is possible that the chemical system is in an apparent chemical equilibrium state when the chemical reactions proceed at extremely slow rates. In such cases, chemical equilibrium might only be achieved after extremely long times (e.g., hundreds to millions of years).

A chemical equilibrium model for the reactions is capable of calculating the final equilibrium state of the chemical system without tracing its intermediate states over time. In general, using chemical equilibrium is a more efficient computational approach for calculating the equilibrium state of the system than using chemical kinetics, which can require many calculation steps until equilibrium is achieved. Many such calculation steps are needed for chemical kinetics because each intermediate state calculated over time must be sufficiently accurate, since future states, including the final equilibrium state, depend on past states. Also note that in a chemical equilibrium model, the initial condition for the chemical state of the system can be given in terms of amounts of components (chemical elements and electrical charge), instead of amounts species, in contrast to a chemical kinetics model. The calculated amounts of each species in the final equilibrium state are constrained, nevertheless, to satisfy the principle of conservation of both chemical element amounts and electrical charge.


---

@todo Finish this section about reactions modeled by a mixed chemical kinetis-equilibrium approach.

---

### Deciding between chemical realism and computational efficiency

One can define a chemical system with many phases, each phase containing one or more species. This does not mean that all phases and their species exist at positive amounts! What it means is that the chemical calculations, equilibrium or kinetics, are capable of deciding if a phase in the chemical system should exist at positive amounts for some given conditions (e.g., temperature, pressure, overall composition).

By selecting as many phases as possible, with the possibilities constrained by the *thermodynamic database* being used, one can increase the confidence level of the estimated chemical states. Note, however, that accurate and realistic estimates depend on many more factors than just the selection of potential phases, such as  the *choice of thermodynamic models for non-ideal phases*. Furthermore, note that adding too many phases and species to the definition of the chemical system can result in *more computationally expensive* chemical calculations. In critical performance applications, such as when combining chemical reactions and fluid flow and species transport modeling, restricting the number of phases and species might be necessary for achieving feasible simulation times. The modeler is responsible to decide to which extent the number of phases and species can be compromised for efficiency reasons at the expense of chemical realism!

### Thermodynamic Assumption for the Phases
In the discussion above, the phases are modeled using *thermodynamics*. A thermodynamic model for a phase presumes that temperature, pressure, and species concentrations are *uniform* within its boundaries. As a result, all other physical and chemical phase properties, such as density, enthalpy, heat capacity, viscosity, which in turn depend on temperature, pressure, and species concentrations, are the same everywhere inside the phase. 

If your modeling problem considers, for example, a fluid with concentration gradients within it, then your problem involves more than just chemical reactions: it also involves *diffusion*, and possibly other transport mechanisms! If these concentration gradients are insignificant for the given time and space scale of your problem, then these phases can be modeled thermodynamically. If not, it might still be possible to subdivide the time and/or space scales of your problem into smaller ones within which a thermodynamic model for the phases are plausible. 

For example, a fluid flowing along a tube or within a porous medium can experience different temperatures, pressures, and possess different species concentrations at different points in space and time. This fluid, as a whole, certainly has no uniform properties within its boundaries. However, by discretizing the larger space in which the fluid flows into several sufficiently small control volumes, or using fine enough discrete grid points, the assumption that the fluid has uniform properties can be established inside those small control volumes or at those grid points.

### Chemical and Thermodynamic Models for the Phases {#chemical-thermo-models}

Specifying the phases and their species is not enough to fully describe a chemical system in the *computational sense*. Every phase in %Reaktoro has two associated models: a *thermodynamic model* and a *chemical model*. These denominations are not standard in the literature, but they are useful in the differentiation of two needed types of models for a phase.

A *thermodynamic model* is a model for the calculation of *standard thermodynamic properties* of the species. Examples include standard Gibbs energies, or standard chemical potentials, standard molar volumes, standard heat capacities, standard enthalpies, and so forth. These models are functions of *temperature* and *pressure* only. Currently, %Reaktoro natively supports only SUPCRT92 databases, which contains parameters for the revised  *Helgeson-Kirkham-Flowers* (HKF) equations of state for the calculation of standard thermodynamic properties for hundreds of aqueous species, at temperatures 0 to 1000 °C and pressures 1 to 5000 bar. The SUPCRT92 databases also contain Maier--Kelly coefficients for the calculation of standard thermodynamic properties of gases and minerals.

A *chemical model* is a model that describes the *non-ideal behavior* of phases. These models not only depend on temperature and pressure, like the thermodynamic models, but also on the amounts of the species in the phase. To be more precise, on the concentrations of these species, which can be calculated from the amounts of the species.

@todo Continue discussion on chemical and thermodynamic models.

### Fundamentals of Chemical Equilibrium
In a chemical equilibrium calculation, the molar amounts of the *chemical species*: @eqc{n=(n_1,\ldots,n_\mathrm{N}),} that correspond to a state of *minimum Gibbs energy* is calculated. In this calculation, the **following constraints** are specified:

* temperature @eq{T} is constrained;
* pressure @eq{P} is constrained; and
* the molar amounts of *elements*, @eq{b=(b_1,\ldots,b_\mathrm{C})}, are constrained.


---

## Chemical Equilibrium Calculations


In this section, we present tutorials for chemical equilibrium calculations. 


### Calculating the equilibrium state of a H@sub{2}O–NaCl–CO@sub{2} system
In this tutorial, we show how to perform an equilibrium calculation in which 1 kg of H@sub{2}O, 0.1 moles of NaCl, and 100 g of CO@sub{2} are mixed at 60 °C and 300 bar. This calculation considers three potential phases: an *aqueous phase*, a *gaseous phase*, and a *mineral phase*. The aqueous phase is defined to represent a saline solution with dissolved CO@sub{2}. The gaseous phase is defined to represent a mixture of gaseous/supercritical carbon dioxide, CO@sub{2}(g), and water vapor, H@sub{2}O(g). The mineral phase is defined to represent halite, NaCl(s), which could precipitate as a result of the aqueous phase becoming saturated with NaCl.

@htmlonly
<a href="tutorial-equilibrium-co2-brine.html"
    <button type="button" class="btn btn-primary">Go to the tutorial »</button>
</a>
@endhtmlonly

<div class="hidden">
@subpage tutorial-equilibrium-co2-brine
</div>

---

### Calculating the equilibrium reaction path of a H@sub{2}O–HCl–CaCO@sub{3} system

@todo Provide more brief details below.

In this tutorial, we show how to perform a sequence of equilibrium calculations that describes the reaction path of a H@sub{2}O--HCl--CaCO@sub{3} system in which HCl is gradually added to the system. 

@htmlonly
<a href="tutorial-equilibriumpath-calcite-hcl.html"
    <button type="button" class="btn btn-primary">Go to the tutorial »</button>
</a>
@endhtmlonly

<div class="hidden">
@subpage tutorial-equilibriumpath-calcite-hcl
</div>

### Chemical Equilibrium Calculations: Inverse Problems

@todo Add here examples in which pH is fixed, and many other constraints.

### Chemical Kinetics Calculations

@todo Add here examples in which chemical kinetics are used.

## Further reading

@todo Improve the organization of these contents.

<!-- - @subpage chemical-equilibrium-calculations -->
<!-- - @subpage chemical-kinetics-calculations -->
- @subpage defining-chemical-systems
- @subpage thermodynamic-databases

