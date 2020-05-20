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
class ChemicalState;
class ChemicalSystem;

//=================================================================================================
//
// Auxiliary Types
//
//=================================================================================================

/// The type of functions implementing equilibrium constraints.
using EquilibriumConstraintFn = Fn<real(const ChemicalProps&)>;

/// The functional equilibrium constraints in a chemical equilibrium calculation.
struct FunctionalConstraints
{
    /// The imposed equilibrium constraint while controling temperature.
    EquilibriumConstraintFn dfdT;

    /// The imposed equilibrium constraint while controling pressure.
    EquilibriumConstraintFn dfdP;

    /// The imposed equilibrium constraints while controlling the titration of substances.
    Pairs<ChemicalFormula, EquilibriumConstraintFn> dfdq;
};

/// The description of a chemical potential constraint in a chemical equilibrium calculation.
struct ChemicalPotentialConstraint
{
    /// The chemical formula of the substance for which the chemical potential is constrained.
    ChemicalFormula formula;

    /// The function of temperature and pressure that evaluates the constrained chemical potential value.
    Fn<real(real,real)> fn;
};

/// The reactivity constraints in a chemical equilibrium calculation.
struct ReactivityConstraints
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
    class Control; class Until; class Attained; class Fix; class Prevent; struct Data;

    /// Construct a EquilibriumConstraints object.
    EquilibriumConstraints(const ChemicalSystem& system);

    /// Destroy this EquilibriumConstraints object.
    ~EquilibriumConstraints();

    /// Return a Control object to initiate the imposition of a general equilibrium constraint.
    auto control() -> Control;

    /// Return a Fix object to initiate the imposition of a chemical potential constraint.
    auto fix() -> Fix;

    /// Return a Prevent object to initiate the imposition of a reactivity constraint.
    auto prevent() -> Prevent;

    /// Return the imposed constraints.
    auto data() const -> const Data&;

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
    /// @param constraints The reference to the underlying functional constraints in EquilibriumConstraints.
    explicit Control(FunctionalConstraints& fconstraints);

    /// Deleted copy constructor.
    Control(const Control&) = delete;

    /// Enable temperature control during the chemical equilibrium calculation.
    auto temperature() -> Until;

    /// Enable pressure control during the chemical equilibrium calculation.
    auto pressure() -> Until;

    /// Enable titration of a substance during the chemical equilibrium calculation.
    auto titrationOf(String titrant) -> Until;

    /// Enable titration of one of the given substances during the chemical equilibrium calculation.
    auto titrationOfEither(String titrant1, String titrant2) -> Until;

private:
    /// The reference to the underlying functional constraints in EquilibriumConstraints.
    FunctionalConstraints& fconstraints;
};

//=================================================================================================
//
// EquilibriumConstraints::Until
//
//=================================================================================================

/// The intermediate auxiliary class used to define the desired equilibrium constraint.
class EquilibriumConstraints::Until
{
public:
    /// Construct an EquilibriumConstraints::Until object.
    /// @param constraintfn The reference to the underlying equilibrium constraint function to be filled in.
    explicit Until(EquilibriumConstraintFn& constraintfn);

    /// Deleted copy constructor.
    Until(const Until&) = delete;

    /// Return an Attained object for setting the underlying equilibrium constraint function.
    auto until() const -> Attained;

    /// Enforce a custom equilibrium constraint at chemical equilibrium.
    auto until(const EquilibriumConstraintFn& customfn) -> void;

private:
    /// The reference to the underlying equilibrium constraint function to be filled in.
    EquilibriumConstraintFn& constraintfn;
};

//=================================================================================================
//
// EquilibriumConstraints::Attained
//
//=================================================================================================

/// The auxiliary class used to define the equilibrium constraint function.
class EquilibriumConstraints::Attained
{
public:
    /// Construct an EquilibriumConstraints::Attained object.
    /// @param constraintfn The reference to the underlying equilibrium constraint function to be filled in.
    explicit Attained(EquilibriumConstraintFn& constraintfn);

    /// Deleted copy constructor.
    Attained(const Attained&) = delete;

    /// Enforce a value for the **volume** of the system at chemical equilibrium.
    /// @param value The constrained volume of the system
    /// @param unit The unit of the constrained volume value (must be convertible to m@sup{3})
    auto volume(real value, String unit) -> void;

    /// Enforce a value for the **internal energy** of the system at chemical equilibrium.
    /// @param value The constrained internal energy of the system
    /// @param unit The unit of the constrained internal energy value (must be convertible to J)
    auto internalEnergy(real value, String unit) -> void;

    /// Enforce a value for the **enthalpy** of the system at chemical equilibrium.
    /// @param value The constrained enthalpy of the system
    /// @param unit The unit of the constrained enthalpy value (must be convertible to J)
    auto enthalpy(real value, String unit) -> void;

    /// Enforce a value for the **Gibbs energy** of the system at chemical equilibrium.
    /// @param value The constrained Gibbs energy of the system
    /// @param unit The unit of the constrained Gibbs energy value (must be convertible to J)
    auto gibbsEnergy(real value, String unit) -> void;

    /// Enforce a value for the **Helmholtz energy** of the system at chemical equilibrium.
    /// @param value The constrained Helmholtz energy of the system
    /// @param unit The unit of the constrained Helmholtz energy value (must be convertible to J)
    auto helmholtzEnergy(real value, String unit) -> void;

    /// Enforce a value for the **entropy** of the system at chemical equilibrium.
    /// @param value The constrained entropy of the system
    /// @param unit The unit of the constrained entropy value (must be convertible to J/K)
    auto entropy(real value, String unit) -> void;

private:
    /// The reference to the underlying equilibrium constraint function to be filled in.
    EquilibriumConstraintFn& constraintfn;
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
    /// Construct a default EquilibriumConstraints::Fix object.
    Fix(const ChemicalSystem& system, ChemicalPotentialConstraint& uconstraint);

    /// Deleted copy constructor.
    Fix(const Fix&) = delete;

    /// Enforce a value for the **chemical potential** of a substance at chemical equilibrium.
    /// @param substance The chemical formula of the substance
    /// @param value The constrained chemical potential value
    /// @param unit The unit for the constrained chemical potential value (must be convertible to J/mol)
    auto chemicalPotential(String substance, real value, String unit) -> void;

    /// Enforce a value for the **ln activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species as found in the database in use
    /// @param value The constrained ln activity value
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto lnActivity(String species, real value) -> void;

    /// Enforce a value for the **lg activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species as found in the database in use
    /// @param value The constrained ln activity value
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto lgActivity(String species, real value) -> void;

    /// Enforce a value for the **activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species as found in the database in use
    /// @param value The constrained activity value
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto activity(String species, real value) -> void;

    /// Enforce a value for the **fugacity** of a gaseous species at chemical equilibrium.
    /// @param species The name of the gaseous species as found in the database in use
    /// @param value The constrained fugacity value
    /// @param unit The unit for the constrained fugacity value (must be convertible to Pa)
    /// @note The gaseous species does not need to be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a gaseous species with given name.
    auto fugacity(String species, real value, String unit) -> void;

    /// Enforce a value for pH at chemical equilibrium.
    /// @param value The constrained value for pH
    auto pH(real value) -> void;

    /// Enforce a value for pMg at chemical equilibrium.
    /// @param value The constrained value for pMg
    auto pMg(real value) -> void;

    /// Enforce a value for pe at chemical equilibrium.
    /// @param value The constrained value for pe
    auto pe(real value) -> void;

    /// Enforce a value for Eh at chemical equilibrium.
    /// @param value The constrained value for Eh
    /// @param unit The unit of the constrained value for Eh (must be convertible to V)
    auto Eh(real value, String unit) -> void;

private:
    /// The chemical system associated with the equilibrium constraints.
    const ChemicalSystem& system;

    /// The substance formula and its associated function that evaluates its constrained chemical potential.
    ChemicalPotentialConstraint& uconstraint;
};

//=================================================================================================
//
// EquilibriumConstraints::Prevent
//
//=================================================================================================

/// The auxiliary class used to define reactivity constraints.
class EquilibriumConstraints::Prevent
{
public:
    /// Construct an EquilibriumConstraints::Prevent object.
    /// @param constraints The reactivity constraints to be filled in by this Prevent object.
    explicit Prevent(const ChemicalSystem& system, ReactivityConstraints& constraints);

    /// Deleted copy constructor.
    Prevent(const Prevent&) = delete;

    /// Prevent the amount of a species from changing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species that should be inert
    auto fromReacting(Index ispecies) -> void;

    /// Prevent a reaction from undergoing any progress during the chemical equilibrium calculation.
    /// @param equation The reaction equation (pairs of species index and stoichiometry) that should be inert
    auto fromReacting(Pairs<Index, double> reaction) -> void;

    /// Prevent a species from reacting or a reaction from undergoing any progress during the chemical equilibrium calculation.
    /// @param what The species name or the reaction equation that should be inert
    auto fromReacting(String what) -> void;

    /// Prevent the amount of a species from increasing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species whose amount cannot increase
    auto fromIncreasing(Index ispecies) -> void;

    /// Prevent the amount of a species from increasing during the chemical equilibrium calculation.
    /// @param species The name of the species whose amount cannot increase
    auto fromIncreasing(String species) -> void;

    /// Prevent the amount of a species from decreasing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species whose amount cannot decrease
    auto fromDecreasing(Index ispecies) -> void;

    /// Prevent the amount of a species from decreasing during the chemical equilibrium calculation.
    /// @param species The name of the species whose amount cannot decrease
    auto fromDecreasing(String species) -> void;

private:
    /// The chemical system associated with the equilibrium constraints.
    const ChemicalSystem& system;

    /// The reactivity constraints to be filled in by this Prevent object.
    ReactivityConstraints& constraints;
};

//=================================================================================================
//
// EquilibriumConstraints::Data
//
//=================================================================================================

/// The auxiliary struct used to store the imposed equilibrium constraints.
struct EquilibriumConstraints::Data
{
    /// The functional equilibrium constraints imposed via method @ref EquilibriumConstraints::control.
    FunctionalConstraints fconstraints;

    /// The description of a chemical potential constraint imposed via method @ref EquilibriumConstraints::fix.
    Vec<ChemicalPotentialConstraint> uconstraints;

    /// The reactivity constraints imposed via method @ref EquilibriumConstraints::prevent.
    ReactivityConstraints rconstraints;
};

} // namespace Reaktoro
