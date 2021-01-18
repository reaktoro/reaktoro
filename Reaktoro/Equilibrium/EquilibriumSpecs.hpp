// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class Params;

/// The details of an equation constraint in a chemical equilibrium calculation.
struct EquilibriumConstraintEquation
{
    /// The unique name of this equation constraint.
    String name;

    /// The function defining the equation to be satisfied at chemical equilibrium.
    /// @param props The chemical properties of the system.
    /// @param params The parameters of the chemical equilibrium calculation.
    Fn<real(const ChemicalProps& props, const Params& params)> fn;
};

/// The details of a chemical potential constraint in a chemical equilibrium calculation.
struct EquilibriumConstraintChemicalPotential
{
    /// The unique name of this chemical potential constraint.
    String name;

    /// The chemical formula of the substance for which the chemical potential is constrained.
    ChemicalFormula substance;

    /// The function that evaluates the constrained chemical potential value.
    /// @param props The chemical properties of the system.
    /// @param params The parameters of the chemical equilibrium calculation.
    Fn<real(const ChemicalProps& props, const Params& params)> fn;
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
/// specs.temperature(); // temperature is an input parameter in the chemical equilibrium problem.
/// specs.pressure();    // pressure is an input parameter in the chemical equilibrium problem.
/// specs.volume();      // volume is an input parameter in the chemical equilibrium problem.
/// specs.openTo("CO2"); // a titrant named [CO2] is introduced, and its amount computed at the end of the equilibrium calculation.
/// ~~~
class EquilibriumSpecs
{
public:
    /// Construct an EquilibriumSpecs object.
    explicit EquilibriumSpecs(const ChemicalSystem& system);

    //=================================================================================================
    //
    // METHODS TO SPECIFY THERMODYNAMIC CONSTRAINTS
    //
    //=================================================================================================

    /// Specify that the **temperature** of the system at chemical equilibrium is given.
    /// This method introduces one input parameter with name `T`. By calling
    /// this method, you are specifying that temperature is known in the
    /// equilibrium calculation. Thus, temperature will not be considered as a
    /// control variable whose value needs to be computed.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesControlVariables.
    auto temperature() -> void;

    /// Specify that the **pressure** of the system at chemical equilibrium is given.
    /// This method introduces one input parameter with name `P`. By calling
    /// this method, you are specifying that pressure is known in the
    /// equilibrium calculation. Thus, pressure will not be considered as a
    /// control variable whose value needs to be computed.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesControlVariables.
    auto pressure() -> void;

    /// Specify that the **volume** of the system at chemical equilibrium is given.
    /// This method introduces one input parameter with name `V`. It also
    /// introduces an equation constraint with name `volume` to enforce a given
    /// volume value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints.
    auto volume() -> void;

    /// Specify that the **internal energy** of the system at chemical equilibrium is given.
    /// This method introduces one input parameter with name `U`. It also
    /// introduces an equation constraint with name `internalEnergy` to enforce a given
    /// internal energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints.
    auto internalEnergy() -> void;

    /// Specify that the **enthalpy** of the system at chemical equilibrium is given.
    /// This method introduces one input parameter with name `H`. It also
    /// introduces an equation constraint with name `enthalpy` to enforce a given
    /// enthalpy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints.
    auto enthalpy() -> void;

    /// Specify that the **Gibbs energy** of the system at chemical equilibrium is given.
    /// This method introduces one input parameter with name `G`. It also
    /// introduces an equation constraint with name `gibbsEnergy` to enforce a given
    /// Gibbs energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints.
    auto gibbsEnergy() -> void;

    /// Specify that the **Helmholtz energy** of the system at chemical equilibrium is given.
    /// This method introduces one input parameter with name `A`. It also
    /// introduces an equation constraint with name `helmholtzEnergy` to enforce a given
    /// Helmholtz energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints.
    auto helmholtzEnergy() -> void;

    /// Specify that the **entropy** of the system at chemical equilibrium is given.
    /// This method introduces one input parameter with name `S`. It also
    /// introduces an equation constraint with name `entropy` to enforce a given
    /// entropy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints.
    auto entropy() -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY CHEMICAL POTENTIAL CONSTRAINTS
    //
    //=================================================================================================

    /// Specify that the **chemical potential** of a substance at chemical equilibrium is given.
    /// This method introduces one input parameter with name `u[substance]`
    /// (e.g., `u[H2O]` if @p substance is `"H2O"`). It also introduces a
    /// chemical potential constraint with same name to enforce a given
    /// chemical potential value for the substance at chemical equilibrium.
    /// This method also indicates that the chemical system is open to @p
    /// substance (e.g., the system is open to mass in/out of H@sub{2}O). Thus
    /// an *implicit titrant* is introduced with name `[substance]` (e.g., `[H2O]`).
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param substance The chemical formula of the substance (e.g., `H2O`, `CO2`, `H+`, `Mg+2`).
    auto chemicalPotential(String substance) -> void;

    /// Specify that the **ln activity** of a species at chemical equilibrium is given.
    /// This method introduces one input parameter with name
    /// `lnActivity[speciesName]` (e.g., `lnActivity[CO2(aq)]` if @p species is
    /// a Species object with name `"CO2(aq)"`). It also introduces a chemical
    /// potential constraint with same name that is equivalent to enforcing a
    /// given value for the natural log activity of the species at chemical
    /// equilibrium. This method also indicates that the chemical system is
    /// open to the underlying substance of the species, not the species
    /// itself. For example, the system is open to mass in/out of CO@sub{2} if
    /// the Species object @p species is `CO2(aq)`, `CO2(g)` or `CO2(l)`). Thus
    /// an *implicit titrant* is introduced with name `[substance]` (e.g., `[CO2]`).
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The chemical species as an Species object.
    auto lnActivity(const Species& species) -> void;

    /// Specify that the **ln activity** of a species at chemical equilibrium is given.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(const Species&).
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto lnActivity(String species) -> void;

    /// Specify that the **lg activity** of a species at chemical equilibrium is given.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(const Species&).
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto lgActivity(String species) -> void;

    /// Specify that the **activity** of a species at chemical equilibrium is given.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(const Species&).
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto activity(String species) -> void;

    /// Specify that the **fugacity** of a gaseous species at chemical equilibrium is given.
    /// This method introduces one input parameter with name `f[gas]` (e.g.,
    /// `f[O2]` if @p gas is `"O2"`). It also introduces a chemical potential
    /// constraint with same name that is equivalent to enforcing a given value
    /// for the fugacity of the gas at chemical equilibrium. This method also
    /// indicates that the chemical system is open to @p gas (e.g., the system is
    /// open to mass in/out of O@sub{2}). Thus an *implicit titrant* is
    /// introduced with name `[substance]` (e.g., `[O2]`).
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param gas The name of the gaseous species as found in the database in use.
    /// @note The gaseous species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a gaseous species with given name.
    auto fugacity(String gas) -> void;

    /// Specify that the pH at chemical equilibrium is given.
    /// This method introduces one input parameter with name `pH`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pH at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of H@sup{+}. Thus an *implicit titrant* is introduced with name `[H+]`.
    /// The code below demonstrate the use of this method and its effect on
    /// the list of parameters, titrants and control variables.
    ///
    /// ~~~{.cpp}
    /// using namespace Reaktoro;
    /// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
    /// specs.enthalpy();                     // introduces parameter `H` and constraint `enthalpy`
    /// specs.pressure();                     // introduces parameter `P`
    /// specs.pH();                           // introduces parameter `pH`, constraint `pH`, and titrant `[H+]`
    /// print(specs.namesParameters());       // H, P, pH
    /// print(specs.namesTitrants());         // [H+]
    /// print(specs.namesControlVariables()); // T, n[H+]
    /// print(specs.namesContraints());       // enthalpy, pH
    /// ~~~
    ///
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @warning An error will be thrown if the system does not contain an aqueous species with formula `H+`.
    auto pH() -> void;

    /// Specify that pMg at chemical equilibrium is given.
    /// This method introduces one input parameter with name `pMg`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pMg at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of Mg@sup{2+}. Thus an *implicit titrant* is introduced with name `[Mg+2]`.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @warning An error will be thrown if the system does not contain an aqueous species with formula `Mg+2`.
    auto pMg() -> void;

    /// Specify that pE at chemical equilibrium is given.
    /// This method introduces one input parameter with name `pE`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pE at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of e@sup{-} (the electron substance). Thus an *implicit titrant*
    /// is introduced with name `[e-]`.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    auto pE() -> void;

    /// Specify that Eh at chemical equilibrium is given.
    /// This method introduces one input parameter with name `Eh`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for Eh at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of e@sup{-} (the electron substance). Thus an *implicit titrant*
    /// is introduced with name `[e-]`.
    /// @see EquilibriumSpecs::namesParameters, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
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
    auto openTo(const ChemicalFormula& substance) -> void;

    //=================================================================================================
    //
    // METHODS TO GET THE NUMBER OF INTRODUCED CONSTRAINTS, PARAMETERS, AND CONTROL VARIABLES
    //
    //=================================================================================================

    /// Return the number of all introduced parameters.
    auto numParameters() const -> Index;

    /// Return the number of all introduced control variables.
    auto numControlVariables() const -> Index;

    /// Return the number of all introduced explicit and implicit titrants.
    auto numTitrants() const -> Index;

    /// Return the number of all introduced explicit titrants.
    auto numTitrantsExplicit() const -> Index;

    /// Return the number of all introduced implicit titrants.
    auto numTitrantsImplicit() const -> Index;

    /// Return the number of all introduced equation and chemical potential constraints.
    auto numConstraints() const -> Index;

    /// Return the number of all introduced constraints of equation type.
    auto numConstraintsEquationType() const -> Index;

    /// Return the number of all introduced constraints of fixed chemical potential type.
    auto numConstraintsChemicalPotentialType() const -> Index;

    //=================================================================================================
    //
    // METHODS TO GET THE NAMES OF INTRODUCED CONSTRAINTS, PARAMETERS, AND CONTROL VARIABLES
    //
    //=================================================================================================

    /// Return the names of all introduced parameters.
    auto namesParameters() const -> Strings;

    /// Return the names of all introduced control variables.
    auto namesControlVariables() const -> Strings;

    /// Return the names of all introduced explicit and implicit titrants.
    auto namesTitrants() const -> Strings;

    /// Return the names of all introduced explicit titrants.
    auto namesTitrantsExplicit() const -> Strings;

    /// Return the names of all introduced implicit titrants.
    auto namesTitrantsImplicit() const -> Strings;

    /// Return the names of all introduced equation and chemical potential constraints.
    auto namesConstraints() const -> Strings;

    /// Return the names of all introduced constraints of equation type.
    auto namesConstraintsEquationType() const -> Strings;

    /// Return the names of all introduced constraints of fixed chemical potential type.
    auto namesConstraintsChemicalPotentialType() const -> Strings;

    //=================================================================================================
    //
    // METHODS TO ADD CONSTRAINTS AND PARAMETERS
    //
    //=================================================================================================

    /// Add a new equation constraint to be satisfied at chemical equilibrium.
    auto addConstraint(const EquilibriumConstraintEquation& constraint) -> void;

    /// Add a new chemical potential constraint to be satisfied at chemical equilibrium.
    auto addConstraint(const EquilibriumConstraintChemicalPotential& constraint) -> void;

    /// Add a new input parameter for the chemical equilibrium problem.
    auto addParameter(String param) -> void;

    //=================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================

    /// Return the chemical system associated with the equilibrium conditions.
    auto system() const -> const ChemicalSystem&;

    /// Return true if temperature is unknown in the chemical equilibrium specifications.
    auto isTemperatureUnknown() const -> bool;

    /// Return true if pressure is unknown in the chemical equilibrium specifications.
    auto isPressureUnknown() const -> bool;

    /// Return the chemical formulas of the explicit and implicit titrant substances.
    auto titrants() const -> Vec<ChemicalFormula>;

    /// Return the chemical formulas of the explicit titrant substances.
    auto titrantsExplicit() const -> Vec<ChemicalFormula>;

    /// Return the chemical formulas of the implicit titrant substances.
    auto titrantsImplicit() const -> Vec<ChemicalFormula>;

    /// Return the equation constraints to be satisfied at chemical equilibrium.
    auto constraintsEquationType() const -> Vec<EquilibriumConstraintEquation> const&;

    /// Return the chemical potential constraints to be satisfied at chemical equilibrium.
    auto constraintsChemicalPotentialType() const -> Vec<EquilibriumConstraintChemicalPotential> const&;

private:
    /// The chemical system associated with the equilibrium conditions.
    ChemicalSystem msystem;

    /// The boolean flag that indicates whether temperature is unknown.
    bool unknownT = true;

    /// The boolean flag that indicates whether pressure is unknown.
    bool unknownP = true;

    /// The input parameters in the chemical equilibrium calculation.
    Strings parameters;

    /// The chemical formulas of the explicit titrants whose amounts are unknown.
    Vec<ChemicalFormula> titrants_explicit;

    /// The chemical formulas of the implicit titrants whose amounts are unknown.
    Vec<ChemicalFormula> titrants_implicit;

    /// The equation constraints to be satisfied at chemical equilibrium.
    Vec<EquilibriumConstraintEquation> econstraints;

    /// The chemical potential constraints to be satisfied at chemical equilibrium.
    Vec<EquilibriumConstraintChemicalPotential> uconstraints;

private:
    /// Throw an error if a given titrant has already been registered explicitly or implicitly.
    auto throwErrorIfTitrantHasBeenRegistered(const ChemicalFormula& substance) const -> void;
};

} // namespace Reaktoro
