// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Params.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;

/// Used to define a *q* control variable in a chemical equilibrium problem.
/// *q* control variables are used to specify the chemical potential of a
/// species at equilibrium. Its introduction in a chemical equilibrium problem
/// implies that the system is open to a titrant with **same** formula of the
/// species whose chemical potential needs to be prescribed. A *q* control
/// variable can also be used to specify the activity of a species. It can also
/// be used to prescribe pH. This can be achieved by converting the given
/// activity value (or pH value) to a chemical potential value of the species.
struct ControlVariableQ
{
    /// The signature of functions that evaluate the prescribed chemical potential of a substance.
    /// @param state The current chemical state and properties of the system in the equilibrium calculation.
    /// @param w The input variables in the chemical equilibrium calculation.
    using ChemicalPotentialFn = Fn<real(ChemicalState const& state, VectorXrConstRef w)>;

    /// The unique name for this *q* control variable (required).
    String name;

    /// The chemical formula of the substance associated to this *q* control variable (required).
    ChemicalFormula substance;

    /// The unique identifier for the chemical potential constraint associated to this *q* control variable (required).
    String id;

    /// The chemical potential function associated to this *q* control variable (required).
    ChemicalPotentialFn fn;
};

/// Used to define a *p* control variable in a chemical equilibrium problem.
/// *p* control variables are more versatile than *q* control variables. They
/// can be used to designate temperature and/or pressure as unknowns. They can
/// be used to represent the amounts of substances (titrants) for which the
/// chemical system is open to. And they can also be used to specify that the
/// chemical potential of a species or a model parameter are unknowns to be
/// resolved along with all others.
struct ControlVariableP
{
    /// The signature of functions that evaluate the chemical potential of a species in terms of a *p* control variable.
    /// @param state The current chemical state and properties of the system in the equilibrium calculation.
    /// @param pk The current value of the corresponding *p* control variable during the equilibrium calculation.
    using ChemicalPotentialFn = Fn<real(ChemicalState const& state, real const& pk)>;

    /// The unique name for this *p* control variable (required).
    String name;

    /// The chemical formula of the substance associated to this *p* control variable (optional).
    ChemicalFormula substance;

    /// The index of the species whose chemical potential is unknown and defined in terms of this *p* control variable (optional).
    Index ispecies = Index(-1);

    /// The function that introduces association between this *p* control variable and the chemical potential of a species (optional).
    ChemicalPotentialFn fn;
};

/// Used to define equation constraints in a chemical equilibrium problem.
struct ConstraintEquation
{
    /// The signature of functions that evaluate the residual of the equation constraint.
    /// @param state The current chemical state and properties of the system in the equilibrium calculation.
    /// @param w The input variables in the chemical equilibrium calculation.
    using ConstraintFn = Fn<real(ChemicalState const& state, VectorXrConstRef w)>;

    /// The unique identifier for this equation constraint.
    String id;

    /// The function defining the equation to be satisfied at chemical equilibrium.
    ConstraintFn fn;
};

/// Used to define reactivity restrictions among species in the chemical
/// equilibrium calculation. This can be used, for example, to impose that a
/// reaction is inert and should not progress during the equilibration process.
struct ReactivityConstraint
{
    /// The unique identifier for this reactivity constraint.
    String id;

    /// The linear equation coefficients in the constraint corresponding to the species amounts variables *n*.
    VectorXd An;

    /// The linear equation coefficients in the constraint corresponding to the introduced control variables *p*.
    VectorXd Ap;
};

/// The class used to define conditions to be satisfied at chemical equilibrium.
///
/// ### Explicit Titrants
///
/// The *explicit titrants* are the titrants introduced with method
/// EquilibriumSpecs::openTo. In the code below, the names of the explicitly
/// introduced titrants are shown:
/// ~~~{.cpp}
/// using namespace Reaktoro;
/// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
/// specs.openTo("H2S"); // a titrant named [H2S] will be introduced.
/// specs.openTo("CO2"); // a titrant named [CO2] will be introduced.
/// ~~~
/// The amounts of these explicit titrants are unknown in the chemical
/// equilibrium problem and computed together with the amounts of species.
///
/// ### Implicit Titrants
///
/// The *implicit titrants* are the titrants introduced with methods:
///
/// * EquilibriumSpecs::chemicalPotential
/// * EquilibriumSpecs::lnActivity
/// * EquilibriumSpecs::lgActivity
/// * EquilibriumSpecs::activity
/// * EquilibriumSpecs::fugacity
/// * EquilibriumSpecs::pH
/// * EquilibriumSpecs::pMg
/// * EquilibriumSpecs::pE
/// * EquilibriumSpecs::Eh
///
/// In the code below, the names of the implicitly introduced titrants are shown:
///
/// ~~~{.cpp}
/// using namespace Reaktoro;
/// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
/// specs.lnActivity("Ca+2"); // a titrant named [Ca+2] will be introduced.
/// specs.fugacity("O2");     // a titrant named [O2] will be introduced.
/// specs.pH();               // a titrant named [H+] will be introduced.
/// specs.pMg();              // a titrant named [Mg+2] will be introduced.
/// specs.pE();               // a titrant named [e-] will be introduced.
/// specs.Eh();               // a titrant named [e-] will be introduced.
/// ~~~
///
/// The amounts of these implicit titrants are unknown in the chemical
/// equilibrium problem and computed together with the amounts of species.
///
/// ### Control Variables
///
/// The *control variables* in a chemical equilibrium problem are unknowns
/// introduced along with equilibrium constraints. These control variables can be:
///
/// * temperature,
/// * pressure, and
/// * amounts of explicit and implicit titrants.
///
/// The number of these control variables depend whether temperature and/or
/// pressure are unknown and if any titrant has been introduced, explicitly
/// or implicitly. For example, if EquilibriumSpecs::temperature is not
/// called, then temperature is unknown and computed in the chemical
/// equilibrium calculation. The same applies for pressure in case
/// EquilibriumSpecs::pressure is not called. Titrants are introduced
/// either explicitly, with method EquilibriumSpecs::openTo, or implicitly,
/// with methods:
///
/// * EquilibriumSpecs::chemicalPotential
/// * EquilibriumSpecs::lnActivity
/// * EquilibriumSpecs::lgActivity
/// * EquilibriumSpecs::activity
/// * EquilibriumSpecs::fugacity
/// * EquilibriumSpecs::pH
/// * EquilibriumSpecs::pMg
/// * EquilibriumSpecs::pE
/// * EquilibriumSpecs::Eh
///
/// The following example formulates a set of specifications for a chemical
/// equilibrium problem in which temperature and the amount of titrant
/// `[CO2]` are introduced *control variables*. Their values are not known
/// in advance; they are computed as part of the chemical equilibrium calculation.
///
/// ~~~{.cpp}
/// using namespace Reaktoro;
/// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
/// specs.temperature(); // temperature is an input variable in the chemical equilibrium problem.
/// specs.pressure();    // pressure is an input variable in the chemical equilibrium problem.
/// specs.volume();      // volume is an input variable in the chemical equilibrium problem.
/// specs.openTo("CO2"); // a titrant named [CO2] is introduced, and its amount computed at the end of the equilibrium calculation.
/// ~~~
class EquilibriumSpecs
{
public:
    /// Construct an EquilibriumSpecs object.
    explicit EquilibriumSpecs(ChemicalSystem const& system);

    //=================================================================================================
    //
    // METHODS TO SPECIFY THERMODYNAMIC CONSTRAINTS
    //
    //=================================================================================================

    /// Specify that the **temperature** of the system at chemical equilibrium is given.
    /// This method introduces one input variable with name `T`. By calling
    /// this method, you are specifying that temperature is known in the
    /// equilibrium calculation. Thus, temperature will not be considered as a
    /// control variable whose value needs to be computed.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesControlVariables.
    auto temperature() -> void;

    /// Specify that the **pressure** of the system at chemical equilibrium is given.
    /// This method introduces one input variable with name `P`. By calling
    /// this method, you are specifying that pressure is known in the
    /// equilibrium calculation. Thus, pressure will not be considered as a
    /// control variable whose value needs to be computed.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesControlVariables.
    auto pressure() -> void;

    /// Specify that the **volume** of the system at chemical equilibrium is given.
    /// This method introduces one input variable with name `V`. It also
    /// introduces an equation constraint with name `volume` to enforce a given
    /// volume value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints.
    auto volume() -> void;

    /// Specify that the **internal energy** of the system at chemical equilibrium is given.
    /// This method introduces one input variable with name `U`. It also
    /// introduces an equation constraint with name `internalEnergy` to enforce a given
    /// internal energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints.
    auto internalEnergy() -> void;

    /// Specify that the **enthalpy** of the system at chemical equilibrium is given.
    /// This method introduces one input variable with name `H`. It also
    /// introduces an equation constraint with name `enthalpy` to enforce a given
    /// enthalpy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints.
    auto enthalpy() -> void;

    /// Specify that the **Gibbs energy** of the system at chemical equilibrium is given.
    /// This method introduces one input variable with name `G`. It also
    /// introduces an equation constraint with name `gibbsEnergy` to enforce a given
    /// Gibbs energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints.
    auto gibbsEnergy() -> void;

    /// Specify that the **Helmholtz energy** of the system at chemical equilibrium is given.
    /// This method introduces one input variable with name `A`. It also
    /// introduces an equation constraint with name `helmholtzEnergy` to enforce a given
    /// Helmholtz energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints.
    auto helmholtzEnergy() -> void;

    /// Specify that the **entropy** of the system at chemical equilibrium is given.
    /// This method introduces one input variable with name `S`. It also
    /// introduces an equation constraint with name `entropy` to enforce a given
    /// entropy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints.
    auto entropy() -> void;

    /// Specify that the **electric charge** at chemical equilibrium is given.
    auto charge() -> void;

    /// Specify that the **amount of an element** at chemical equilibrium is given.
    /// @param element The name or index of the element in the chemical system.
    auto elementAmount(StringOrIndex const& element) -> void;

    /// Specify that the **amount of an element in a phase** at chemical equilibrium is given.
    /// @param element The name or index of the element in the chemical system.
    /// @param phase The name or index of the phase in the chemical system.
    auto elementAmountInPhase(StringOrIndex const& element, StringOrIndex const& phase) -> void;

    /// Specify that the **mass of an element** at chemical equilibrium is given.
    /// @param element The name or index of the element in the chemical system.
    auto elementMass(StringOrIndex const& element) -> void;

    /// Specify that the **mass of an element in a phase** at chemical equilibrium is given.
    /// @param element The name or index of the element in the chemical system.
    /// @param phase The name or index of the phase in the chemical system.
    auto elementMassInPhase(StringOrIndex const& element, StringOrIndex const& phase) -> void;

    /// Specify that the **amount of a phase** at chemical equilibrium is given.
    /// @param phase The name or index of the phase in the chemical system.
    auto phaseAmount(StringOrIndex const& phase) -> void;

    /// Specify that the **mass of a phase** at chemical equilibrium is given.
    /// @param phase The name or index of the phase in the chemical system.
    auto phaseMass(StringOrIndex const& phase) -> void;

    /// Specify that the **volume of a phase** at chemical equilibrium is given.
    /// @param phase The name or index of the phase in the chemical system.
    auto phaseVolume(StringOrIndex const& phase) -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY CHEMICAL POTENTIAL CONSTRAINTS
    //
    //=================================================================================================

    /// Specify that the **chemical potential** of a substance at chemical equilibrium is given.
    /// This method introduces one input variable with name `u[substance]`
    /// (e.g., `u[H2O]` if @p substance is `"H2O"`). It also introduces a
    /// chemical potential constraint with same name to enforce a given
    /// chemical potential value for the substance at chemical equilibrium.
    /// This method also indicates that the chemical system is open to @p
    /// substance (e.g., the system is open to mass in/out of H@sub{2}O). Thus
    /// an *implicit titrant* is introduced with name `[substance]` (e.g., `[H2O]`).
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param substance The chemical formula of the substance (e.g., `H2O`, `CO2`, `H+`, `Mg+2`).
    auto chemicalPotential(String substance) -> void;

    /// Specify that the **ln activity** of a species at chemical equilibrium is given.
    /// This method introduces one input variable with name
    /// `lnActivity[speciesName]` (e.g., `lnActivity[CO2(aq)]` if @p species is
    /// a Species object with name `"CO2(aq)"`). It also introduces a chemical
    /// potential constraint with same name that is equivalent to enforcing a
    /// given value for the natural log activity of the species at chemical
    /// equilibrium. This method also indicates that the chemical system is
    /// open to the underlying substance of the species, not the species
    /// itself. For example, the system is open to mass in/out of CO@sub{2} if
    /// the Species object @p species is `CO2(aq)`, `CO2(g)` or `CO2(l)`). Thus
    /// an *implicit titrant* is introduced with name `[substance]` (e.g., `[CO2]`).
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The chemical species as an Species object.
    auto lnActivity(Species const& species) -> void;

    /// Specify that the **ln activity** of a species at chemical equilibrium is given.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(Species const&).
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto lnActivity(String species) -> void;

    /// Specify that the **lg activity** of a species at chemical equilibrium is given.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(Species const&).
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto lgActivity(String species) -> void;

    /// Specify that the **activity** of a species at chemical equilibrium is given.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(Species const&).
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto activity(String species) -> void;

    /// Specify that the **fugacity** of a gaseous species at chemical equilibrium is given.
    /// This method introduces one input variable with name `f[gas]` (e.g.,
    /// `f[O2]` if @p gas is `"O2"`). It also introduces a chemical potential
    /// constraint with same name that is equivalent to enforcing a given value
    /// for the fugacity of the gas at chemical equilibrium. This method also
    /// indicates that the chemical system is open to @p gas (e.g., the system is
    /// open to mass in/out of O@sub{2}). Thus an *implicit titrant* is
    /// introduced with name `[substance]` (e.g., `[O2]`).
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param gas The name of the gaseous species as found in the database in use.
    /// @note The gaseous species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a gaseous species with given name.
    auto fugacity(String gas) -> void;

    /// Specify that the pH at chemical equilibrium is given.
    /// This method introduces one input variable with name `pH`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pH at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of H@sup{+}. Thus an *implicit titrant* is introduced with name `[H+]`.
    /// The code below demonstrate the use of this method and its effect on
    /// the list of input variables, titrants and control variables.
    ///
    /// ~~~{.cpp}
    /// using namespace Reaktoro;
    /// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
    /// specs.enthalpy();                     // introduces input variable `H` and constraint `enthalpy`
    /// specs.pressure();                     // introduces input variable `P`
    /// specs.pH();                           // introduces input variable `pH`, constraint `pH`, and titrant `[H+]`
    /// print(specs.namesParams());           // H, P, pH
    /// print(specs.namesTitrants());         // [H+]
    /// print(specs.namesControlVariables()); // T, n[H+]
    /// print(specs.namesContraints());       // enthalpy, pH
    /// ~~~
    ///
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @warning An error will be thrown if the system does not contain an aqueous species with formula `H+`.
    auto pH() -> void;

    /// Specify that pMg at chemical equilibrium is given.
    /// This method introduces one input variable with name `pMg`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pMg at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of Mg@sup{2+}. Thus an *implicit titrant* is introduced with name `[Mg+2]`.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @warning An error will be thrown if the system does not contain an aqueous species with formula `Mg+2`.
    auto pMg() -> void;

    /// Specify that pE at chemical equilibrium is given.
    /// This method introduces one input variable with name `pE`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pE at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of e@sup{-} (the electron substance). Thus an *implicit titrant*
    /// is introduced with name `[e-]`.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    auto pE() -> void;

    /// Specify that Eh at chemical equilibrium is given.
    /// This method introduces one input variable with name `Eh`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for Eh at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of e@sup{-} (the electron substance). Thus an *implicit titrant*
    /// is introduced with name `[e-]`.
    /// @see EquilibriumSpecs::namesParams, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    auto Eh() -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY HOW THE CHEMICAL SYSTEM IS OPEN
    //
    //=================================================================================================

    /// Specify that the chemical system is open to a substance.
    ///
    /// Use this method to specify that the system is titrated with an unknown
    /// amount of a substance to be able to attain chemical equilibrium with
    /// given conditions. Its use introduces an *explicit titrant* with name
    /// `[substance]` (e.g., `H2S` is @p substance is `"H2S"`). The amount of
    /// this titrant is an unknown control variable which is computed in the
    /// chemical equilibrium calculation.
    ///
    /// The code below demonstrate the use of this method and its effect on
    /// the list of titrants and control variables.
    ///
    /// ~~~{.cpp}
    /// using namespace Reaktoro;
    /// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
    /// specs.volume();
    /// specs.pressure();
    /// specs.openTo("H2S");
    /// specs.openTo("CO2");
    /// print(specs.namesTitrants()); // [H2S], [CO2]
    /// print(specs.namesControlVariables()); // T, [H2S], [CO2]
    /// ~~~
    ///
    /// @note The code above is for demonstration purposes only. The given
    /// specifications do not produce a valid chemical equilibrium problem!
    ///
    /// @see EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param substance The chemical formula of the substance.
    auto openTo(ChemicalFormula const& substance) -> void;

    /// Specify that the chemical system is open to a given chemical state. // TODO: Implement EquilibriumSpecs::auto openTo(ChemicalState const& state).
    // auto openTo(ChemicalState const& state) -> void;

    /// Specify that the chemical system is open to a given material. // TODO: Implement EquilibriumSpecs::auto openTo(Material const& material).
    // auto openTo(Material const& material) -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY ADDITIONAL UNKNOWNS
    //
    //=================================================================================================

    /// Specify that the chemical system is open to a titrant substance and its amount is unknown.
    auto addUnknownTitrantAmount(ChemicalFormula const& substance) -> void;

    /// Specify that the chemical potential of a species is unknown at equilibrium and must be computed.
    auto addUnknownChemicalPotential(String const& species) -> void;

    /// Specify that the standard chemical potential of a species is unknown at equilibrium and must be computed.
    auto addUnknownStandardChemicalPotential(String const& species) -> void;

    /// Specify that the activity of a species is unknown at equilibrium and must be computed.
    auto addUnknownActivity(String const& species) -> void;

    /// Specify that the activity coefficient of a species is unknown at equilibrium and must be computed.
    auto addUnknownActivityCoefficient(String const& species) -> void;

    //=================================================================================================
    //
    // METHODS TO GET THE NUMBER OF INTRODUCED CONSTRAINTS, INPUT VARIABLES, AND CONTROL VARIABLES
    //
    //=================================================================================================

    /// Return the number of introduced input variables.
    auto numInputs() const -> Index;

    /// Return the number of model parameters among the introduced input variables.
    auto numParams() const -> Index;

    /// Return the number of all introduced control variables.
    auto numControlVariables() const -> Index;

    /// Return the number of introduced *p* control variables.
    auto numControlVariablesP() const -> Index;

    /// Return the number of introduced *q* control variables.
    auto numControlVariablesQ() const -> Index;

    /// Return the number of all introduced explicit and implicit titrants.
    auto numTitrants() const -> Index;

    /// Return the number of all introduced explicit titrants.
    auto numTitrantsExplicit() const -> Index;

    /// Return the number of all introduced implicit titrants.
    auto numTitrantsImplicit() const -> Index;

    /// Return the number of all introduced equation and chemical potential constraints.
    auto numConstraints() const -> Index;

    //=================================================================================================
    //
    // METHODS TO GET THE NAMES OF INTRODUCED CONSTRAINTS, INPUT VARIABLES, AND CONTROL VARIABLES
    //
    //=================================================================================================

    /// Return the names of the introduced input variables.
    auto namesInputs() const -> Strings;

    /// Return the names of the model parameters among the input variables.
    auto namesParams() const -> Strings;

    /// Return the names of all introduced control variables.
    auto namesControlVariables() const -> Strings;

    /// Return the names of introduced *p* control variables.
    auto namesControlVariablesP() const -> Strings;

    /// Return the names of introduced *q* control variables.
    auto namesControlVariablesQ() const -> Strings;

    /// Return the names of all introduced explicit and implicit titrants.
    auto namesTitrants() const -> Strings;

    /// Return the names of all introduced explicit titrants.
    auto namesTitrantsExplicit() const -> Strings;

    /// Return the names of all introduced implicit titrants.
    auto namesTitrantsImplicit() const -> Strings;

    /// Return the names of all introduced equation and chemical potential constraints.
    auto namesConstraints() const -> Strings;

    //=================================================================================================
    //
    // METHODS TO ADD CONTROL VARIABLES, CONSTRAINTS, AND INPUT VARIABLES
    //
    //=================================================================================================

    /// Add a *q* control variable in the specification of the chemical equilibrium problem.
    auto addControlVariableQ(ControlVariableQ const& qvar) -> void;

    /// Add a *p* control variable in the specification of the chemical equilibrium problem.
    auto addControlVariableP(ControlVariableP const& pvar) -> void;

    /// Add an equation constraint to be satisfied at chemical equilibrium.
    auto addConstraint(ConstraintEquation const& constraint) -> void;

    /// Add a reactivity constraint to be satisfied at chemical equilibrium.
    auto addReactivityConstraint(ReactivityConstraint const& constraint) -> void;

    /// Add a new input variable for the chemical equilibrium problem with name @p var.
    auto addInput(String const& var) -> Index;

    /// Add model parameter @p param as a new input variable for the chemical equilibrium problem.
    auto addInput(Param const& param) -> Index;

    //=================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================

    /// Return the chemical system associated with the equilibrium conditions.
    auto system() const -> ChemicalSystem const&;

    /// Return the input variables in the chemical equilibrium specifications.
    auto inputs() const -> Strings const&;

    /// Return the model parameters among the input variables.
    auto params() const -> const Vec<Param>&;

    /// Return the indices of the model parameters among the input variables.
    auto indicesParams() const -> const Vec<Index>&; // TODO: Rename to indicesInputParams because there should be another method called indicesUnknownParams.

    /// Return true if temperature is unknown in the chemical equilibrium specifications.
    auto isTemperatureUnknown() const -> bool;

    /// Return true if pressure is unknown in the chemical equilibrium specifications.
    auto isPressureUnknown() const -> bool;

    /// The index of temperature among the *p* control variables or `Index(-1)` if it is a given input.
    auto indexControlVariableTemperature() const -> Index;

    /// The index of pressure among the *p* control variables or `Index(-1)` if it is a given input.
    auto indexControlVariablePressure() const -> Index;

    /// Return the *q* control variables in the chemical equilibrium specifications.
    auto controlVariablesQ() const -> const Vec<ControlVariableQ>&;

    /// Return the *q* control variables in the chemical equilibrium specifications.
    auto controlVariablesP() const -> const Vec<ControlVariableP>&;

    /// Return the chemical formulas of the explicit and implicit titrant substances.
    auto titrants() const -> Vec<ChemicalFormula>;

    /// Return the chemical formulas of the explicit titrant substances.
    auto titrantsExplicit() const -> Vec<ChemicalFormula>;

    /// Return the chemical formulas of the implicit titrant substances.
    auto titrantsImplicit() const -> Vec<ChemicalFormula>;

    /// Return the equation constraints to be satisfied at chemical equilibrium.
    auto equationConstraints() const -> Vec<ConstraintEquation> const&;

    /// Return the introduced reactivity constraints to be satisfied during the equilibrium calculation.
    auto reactivityConstraints() const -> Vec<ReactivityConstraint> const&;

private:
    /// The chemical system associated with the equilibrium conditions.
    ChemicalSystem m_system;

    /// The names of the input variables in the chemical equilibrium calculation.
    Strings m_inputs;

    /// The model parameters among the input variables in the chemical equilibrium calculation.
    Vec<Param> m_params;

    /// The indices of the model parameters among the input variables in the chemical equilibrium calculation.
    Vec<Index> m_params_idxs;

    /// The boolean flag that indicates whether temperature is unknown.
    bool unknownT = true;

    /// The boolean flag that indicates whether pressure is unknown.
    bool unknownP = true;

    /// The *q* control variables in the chemical equilibrium problem.
    Vec<ControlVariableQ> qvars;

    /// The *p* control variables in the chemical equilibrium problem.
    Vec<ControlVariableP> pvars;

    /// The equation constraints to be satisfied at chemical equilibrium.
    Vec<ConstraintEquation> econstraints;

    /// The reactivity constraints to be satisfied at chemical equilibrium.
    Vec<ReactivityConstraint> rconstraints;

    // ----- AUXILIARY DATA MEMBERS ----- //

    /// The chemical formulas of the explicit titrants whose amounts are unknown.
    Vec<ChemicalFormula> titrants_explicit;

    /// The chemical formulas of the implicit titrants whose amounts are unknown.
    Vec<ChemicalFormula> titrants_implicit;

    /// The names of the species whose chemical potentials are unknown.
    Strings species_with_unknown_chemical_potentials;

private:
    /// Throw an error if a given titrant has already been registered explicitly or implicitly.
    auto throwErrorIfTitrantHasBeenRegistered(ChemicalFormula const& substance) const -> void;
};

} // namespace Reaktoro
