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

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalSystem;
class Species;

//=================================================================================================
//
// EquilibriumConstraints
//
//=================================================================================================

/// The class used to define equilibrium constraints.
class EquilibriumConstraints
{
public:
    // Forward declarations
    class Control; class Until; class Fix; class Preserve; class Prevent; struct Data;

    /// Construct an EquilibriumConstraints object.
    EquilibriumConstraints(const ChemicalSystem& system);

    /// Construct a copy of an EquilibriumConstraints object.
    EquilibriumConstraints(const EquilibriumConstraints& other);

    /// Destroy this EquilibriumConstraints object.
    ~EquilibriumConstraints();

    /// Assign a copy of an EquilibriumConstraints object to this.
    auto operator=(EquilibriumConstraints other) -> EquilibriumConstraints&;

    /// Return a Control object to initiate the introduction of **control variables**.
    auto control() -> Control;

    /// Return a Until object to initiate the imposition of an **equation constraint**.
    auto until() -> Until;

    /// Return a Preserve object to initiate the imposition of a **property preservation constraint**.
    auto preserve() -> Preserve;

    /// Return a Fix object to initiate the imposition of a **chemical potential constraint**.
    auto fix() -> Fix;

    /// Return a Prevent object to initiate the imposition of a **reactivity restriction**.
    auto prevent() -> Prevent;

    /// Return the chemical system associated with the equilibrium constraints.
    auto system() const -> const ChemicalSystem&;

    /// Return the data with details of the imposed constraints.
    auto data() const -> const Data&;

    /// Lock this object to prevent certain changes.
    /// Once locked, an EquilibriumConstraints object will not allow:
    ///   1. the introduction of new control variables;
    ///   2. the imposition of new functional and chemical potential constraints;
    ///   3. the introduction of new inert reactions.
    /// Thus, calling methods @ref control and @ref preserve will result in an
    /// error. Method @ref until can still be called, but only to update an
    /// existing constraint, not to introduce a new one. Method @ref prevent
    /// can be called, as long as it does not introduce a new inert reaction
    /// using method `fromReacting`.
    auto lock() -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

//=================================================================================================
//
// EquilibriumConstraints::Control
//
//=================================================================================================

/// The auxiliary class used to define the controlled variable in an equilibrium constraint.
class EquilibriumConstraints::Control
{
public:
    /// Construct an EquilibriumConstraints::Control object.
    /// @param data The underlying data of the EquilibriumConstraints object.
    explicit Control(EquilibriumConstraints::Data& data);

    /// Enable temperature control during the chemical equilibrium calculation.
    auto temperature() -> Control&;

    /// Enable pressure control during the chemical equilibrium calculation.
    auto pressure() -> Control&;

    /// Enable titration control of a substance during the chemical equilibrium calculation.
    auto titrationOf(String titrant) -> Control&;

    /// Enable titration control of one of the given substances during the chemical equilibrium calculation.
    auto titrationOfEither(String titrant1, String titrant2) -> Control&;

private:
    /// The underlying data of the EquilibriumConstraints object.
    EquilibriumConstraints::Data& data;

    /// Deleted copy constructor.
    Control(const Control&) = delete;
};

//=================================================================================================
//
// EquilibriumConstraints::Until
//
//=================================================================================================

/// The arguments of functions defining equation constraints to be satisfied at chemical equilibrium.
struct EquilibriumEquationArgs
{
    /// The current chemical properties of the chemical system during the calculation.
    const ChemicalProps& props;

    /// The current amounts of the controlled titrants during the calculation.
    VectorXrConstRef q;

    /// The chemical formulas of the controlled titrants.
    const Vec<ChemicalFormula>& titrants;

    /// Return the amount of a controlled titrant with given formula.
    auto titrantAmount(const String& formula) const -> real;
};

/// The type of functions defining equation constraints to be satisfied at chemical equilibrium.
using EquilibriumEquationFn = Fn<real(EquilibriumEquationArgs)>;

/// The auxiliary class used to impose equations constraints at chemical equilibrium.
class EquilibriumConstraints::Until
{
public:
    /// Construct an EquilibriumConstraints::Until object.
    /// @param data The underlying data of the EquilibriumConstraints object.
    explicit Until(EquilibriumConstraints::Data& data);

    /// Enforce a value for the **volume** of the system at chemical equilibrium.
    /// @param value The constrained volume of the system
    /// @param unit The unit of the constrained volume value (must be convertible to m@sup{3})
    auto volume(real value, String unit) -> Until&;

    /// Enforce a value for the **internal energy** of the system at chemical equilibrium.
    /// @param value The constrained internal energy of the system
    /// @param unit The unit of the constrained internal energy value (must be convertible to J)
    auto internalEnergy(real value, String unit) -> Until&;

    /// Enforce a value for the **enthalpy** of the system at chemical equilibrium.
    /// @param value The constrained enthalpy of the system
    /// @param unit The unit of the constrained enthalpy value (must be convertible to J)
    auto enthalpy(real value, String unit) -> Until&;

    /// Enforce a value for the **Gibbs energy** of the system at chemical equilibrium.
    /// @param value The constrained Gibbs energy of the system
    /// @param unit The unit of the constrained Gibbs energy value (must be convertible to J)
    auto gibbsEnergy(real value, String unit) -> Until&;

    /// Enforce a value for the **Helmholtz energy** of the system at chemical equilibrium.
    /// @param value The constrained Helmholtz energy of the system
    /// @param unit The unit of the constrained Helmholtz energy value (must be convertible to J)
    auto helmholtzEnergy(real value, String unit) -> Until&;

    /// Enforce a value for the **entropy** of the system at chemical equilibrium.
    /// @param value The constrained entropy of the system
    /// @param unit The unit of the constrained entropy value (must be convertible to J/K)
    auto entropy(real value, String unit) -> Until&;

    /// Enforce a custom equation constraint at chemical equilibrium.
    /// @param id The unique identifier string of this custom equation constraint
    /// @param fn The function implemeting the custom constraint equation
    auto custom(const String& id, const EquilibriumEquationFn& fn) -> Until&;

private:
    /// The underlying data of the EquilibriumConstraints object.
    EquilibriumConstraints::Data& data;

    /// Deleted copy constructor.
    Until(const Until&) = delete;
};

//=================================================================================================
//
// EquilibriumConstraints::Preserve
//
//=================================================================================================

/// The type of functions that compute or retrieve a chemical property.
using ChemicalPropertyFn = Fn<real(const ChemicalProps&)>;

/// The auxiliary class used to impose a property preservation constraint.
class EquilibriumConstraints::Preserve
{
public:
    /// Construct an EquilibriumConstraints::Preserve object.
    /// @param data The underlying data of the EquilibriumConstraints object.
    explicit Preserve(EquilibriumConstraints::Data& data);

    /// Preserve the **volume** of the system during the chemical equilibrium calculation.
    auto volume() -> Preserve&;

    /// Preserve the **internal energy** of the system during the chemical equilibrium calculation.
    auto internalEnergy() -> Preserve&;

    /// Preserve the **enthalpy** of the system during the chemical equilibrium calculation.
    auto enthalpy() -> Preserve&;

    /// Preserve the **Gibbs energy** of the system during the chemical equilibrium calculation.
    auto gibbsEnergy() -> Preserve&;

    /// Preserve the **Helmholtz energy** of the system during the chemical equilibrium calculation.
    auto helmholtzEnergy() -> Preserve&;

    /// Preserve the **entropy** of the system during the chemical equilibrium calculation.
    auto entropy() -> Preserve&;

    /// Preserve a **custom property** of the system during the chemical equilibrium calculation.
    /// @param id The unique identifier string of this custom property preservation constraint
    /// @param fn The function that computes the custom chemical property.
    auto custom(const String& id, const ChemicalPropertyFn& fn) -> Preserve&;

private:
    /// The underlying data of the EquilibriumConstraints object.
    EquilibriumConstraints::Data& data;

    /// Deleted copy constructor.
    Preserve(const Preserve&) = delete;
};

//=================================================================================================
//
// EquilibriumConstraints::Fix
//
//=================================================================================================

/// The class used to define constraints on chemical potentials of substances.
class EquilibriumConstraints::Fix
{
public:
    /// Construct an EquilibriumConstraints::Fix object.
    /// @param system The underlying chemical system of the EquilibriumConstraints object.
    /// @param data The underlying data of the EquilibriumConstraints object.
    Fix(const ChemicalSystem& system, EquilibriumConstraints::Data& data);

    /// Deleted copy constructor.
    Fix(const Fix&) = delete;

    /// Enforce a value for the **chemical potential** of a substance at chemical equilibrium.
    /// @param substance The chemical formula of the substance
    /// @param fn The function of temperature and pressure that computes the constrained chemical potential value (in J/mol)
    auto chemicalPotential(const ChemicalFormula& substance, const Fn<real(real,real)>& fn) -> Fix&;

    /// Enforce a value for the **chemical potential** of a substance at chemical equilibrium.
    /// @param substance The chemical formula of the substance
    /// @param value The constrained chemical potential value
    /// @param unit The unit for the constrained chemical potential value (must be convertible to J/mol)
    auto chemicalPotential(String substance, real value, String unit) -> Fix&;

    /// Enforce a value for the **ln activity** of a species at chemical equilibrium.
    /// @param species The chemical species
    /// @param value The constrained ln activity value
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto lnActivity(const Species& species, real value) -> Fix&;

    /// Enforce a value for the **ln activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species as found in the database in use
    /// @param value The constrained ln activity value
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto lnActivity(String species, real value) -> Fix&;

    /// Enforce a value for the **lg activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species as found in the database in use
    /// @param value The constrained ln activity value
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto lgActivity(String species, real value) -> Fix&;

    /// Enforce a value for the **activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species as found in the database in use
    /// @param value The constrained activity value
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto activity(String species, real value) -> Fix&;

    /// Enforce a value for the **fugacity** of a gaseous species at chemical equilibrium.
    /// @param species The name of the gaseous species as found in the database in use
    /// @param value The constrained fugacity value
    /// @param unit The unit for the constrained fugacity value (must be convertible to Pa)
    /// @note The gaseous species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a gaseous species with given name.
    auto fugacity(String species, real value, String unit) -> Fix&;

    /// Enforce a value for pH at chemical equilibrium.
    /// @param value The constrained value for pH
    auto pH(real value) -> Fix&;

    /// Enforce a value for pMg at chemical equilibrium.
    /// @param value The constrained value for pMg
    auto pMg(real value) -> Fix&;

    /// Enforce a value for pe at chemical equilibrium.
    /// @param value The constrained value for pe
    auto pe(real value) -> Fix&;

    /// Enforce a value for Eh at chemical equilibrium.
    /// @param value The constrained value for Eh
    /// @param unit The unit of the constrained value for Eh (must be convertible to V)
    auto Eh(real value, String unit) -> Fix&;

private:
    /// The underlying chemical system of the EquilibriumConstraints object.
    const ChemicalSystem& system;

    /// The underlying data of the EquilibriumConstraints object.
    EquilibriumConstraints::Data& data;
};

//=================================================================================================
//
// EquilibriumConstraints::Prevent
//
//=================================================================================================

/// The auxiliary class used to define reactivity restrictions.
class EquilibriumConstraints::Prevent
{
public:
    /// Construct an EquilibriumConstraints::Prevent object.
    /// @param system The underlying chemical system of the EquilibriumConstraints object.
    /// @param data The underlying data of the EquilibriumConstraints object.
    Prevent(const ChemicalSystem& system, EquilibriumConstraints::Data& data);

    /// Deleted copy constructor.
    Prevent(const Prevent&) = delete;

    /// Prevent the amount of a species from changing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species that should be inert
    auto fromReacting(Index ispecies) -> Prevent&;

    /// Prevent a reaction from undergoing any progress during the chemical equilibrium calculation.
    /// @param equation The reaction equation (pairs of species index and stoichiometry) that should be inert
    auto fromReacting(Pairs<Index, double> reaction) -> Prevent&;

    /// Prevent a species from reacting or a reaction from undergoing any progress during the chemical equilibrium calculation.
    /// @param what The species name or the reaction equation that should be inert
    auto fromReacting(String what) -> Prevent&;

    /// Prevent the amount of a species from increasing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species whose amount cannot increase
    auto fromIncreasing(Index ispecies) -> Prevent&;

    /// Prevent the amount of a species from increasing during the chemical equilibrium calculation.
    /// @param species The name of the species whose amount cannot increase
    auto fromIncreasing(String species) -> Prevent&;

    /// Prevent the amount of a species from decreasing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species whose amount cannot decrease
    auto fromDecreasing(Index ispecies) -> Prevent&;

    /// Prevent the amount of a species from decreasing during the chemical equilibrium calculation.
    /// @param species The name of the species whose amount cannot decrease
    auto fromDecreasing(String species) -> Prevent&;

private:
    /// The underlying chemical system of the EquilibriumConstraints object.
    const ChemicalSystem& system;

    /// The underlying data of the EquilibriumConstraints object.
    EquilibriumConstraints::Data& data;
};

//=================================================================================================
//
// EquilibriumConstraints::Data
//
//=================================================================================================

/// The auxiliary struct used to store the imposed equilibrium constraints.
struct EquilibriumConstraints::Data
{
    /// Auxiliary struct containing details of the control variables.
    struct Controls
    {
        /// The boolean flag that indicates whether temperature is controlled.
        bool T = false;

        /// The boolean flag that indicates whether pressure is controlled.
        bool P = false;

        /// The chemical formulas of the titrants whose amounts are controlled.
        Vec<ChemicalFormula> titrants;

        /// Return the number of introduced control variables.
        auto size() const -> Index { return titrants.size() + T + P; }
    };

    /// The details of an equation constraint in a chemical equilibrium calculation.
    struct EquationConstraint
    {
        /// The unique identifier string of this equation constraint.
        String id;

        /// The function defining the equation constraint in residual form.
        EquilibriumEquationFn fn;
    };

    /// The details of a property preservation constraint in a chemical equilibrium calculation.
    struct PropertyPreservationConstraint
    {
        /// The unique identifier string of this property preservation constraint.
        String id;

        /// The function that computes or retrieves the preserved chemical property.
        ChemicalPropertyFn fn;
    };

    /// The details of a chemical potential constraint in a chemical equilibrium calculation.
    struct ChemicalPotentialConstraint
    {
        /// The chemical formula of the substance for which the chemical potential is constrained.
        ChemicalFormula formula;

        /// The function of temperature and pressure that evaluates the constrained chemical potential value.
        Fn<real(real,real)> fn;
    };

    /// The reactivity restrictions in a chemical equilibrium calculation.
    struct ReactivityRestrictions
    {
        /// The indices of the species whose amounts cannot change.
        Set<Index> species_cannot_react;

        /// The indices of the species whose amounts cannot increase.
        Set<Index> species_cannot_increase;

        /// The indices of the species whose amounts cannot decrease.
        Set<Index> species_cannot_decrease;

        /// The inert reactions as pairs of species index and its stoichiometric coefficient.
        Vec<Pairs<Index, double>> reactions_cannot_react;
    };

    /// The introduced control variables via method @ref EquilibriumConstraints::control.
    Controls controls;

    /// The equation constraints imposed via method @ref EquilibriumConstraints::until.
    Vec<EquationConstraint> econstraints;

    /// The property preservation constraints imposed via method @ref EquilibriumConstraints::preserve.
    Vec<PropertyPreservationConstraint> pconstraints;

    /// The chemical potential constraints imposed via method @ref EquilibriumConstraints::fix.
    Vec<ChemicalPotentialConstraint> uconstraints;

    /// The reactivity restrictions imposed via method @ref EquilibriumConstraints::prevent.
    ReactivityRestrictions restrictions;

    /// The flag that indicates whether the EquilibriumConstraints object is locked.
    bool locked = false;
};

} // namespace Reaktoro
