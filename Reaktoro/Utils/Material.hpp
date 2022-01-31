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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>

namespace Reaktoro {

/// A type used to represent a material composed of one or more substances.
class Material
{
public:
    /// Construct a Material object.
    Material(ChemicalSystem const& system);

    /// Add a given amount of an existing species in the chemical system to the material.
    /// @param species The name or index of the species in the chemical system.
    /// @param amount The amount of the species (in mol)
    /// @warning An error is thrown if `species` is not the name or index of an existing species in the chemical system.
    auto addSpeciesAmount(StringOrIndex const& species, double amount) -> void;

    /// Add a given amount of an existing species in the chemical system to the material.
    /// @param species The name or index of the species in the chemical system.
    /// @param amount The amount of the species
    /// @param unit The amount unit
    /// @warning An error is thrown if `species` is not the name or index of an existing species in the chemical system.
    auto addSpeciesAmount(StringOrIndex const& species, double amount, Chars unit) -> void;

    /// Add a given mass of an existing species in the chemical system to the material.
    /// @param species The name or index of the species in the chemical system.
    /// @param mass The mass of the species
    /// @param unit The mass unit
    /// @warning An error is thrown if `species` is not the name or index of an existing species in the chemical system.
    auto addSpeciesMass(StringOrIndex const& species, double mass, Chars unit) -> void;

    /// Add a given amount of a substance to the material.
    /// @param substance The chemical formula of the substance.
    /// @param amount The amount of the substance (in mol)
    auto addSubstanceAmount(ChemicalFormula const& substance, double amount) -> void;

    /// Add a given amount of a substance to the material.
    /// @param substance The chemical formula of the substance.
    /// @param amount The amount of the substance
    /// @param unit The amount unit
    auto addSubstanceAmount(ChemicalFormula const& substance, double amount, Chars unit) -> void;

    /// Add a given mass of a substance to the material.
    /// @param substance The chemical formula of the substance.
    /// @param mass The mass of the substance
    /// @param unit The mass unit
    auto addSubstanceMass(ChemicalFormula const& substance, double mass, Chars unit) -> void;

    /// Add a given amount of another material in this material.
    /// @param material The other material composing this material.
    /// @param amount The amount of the substance (in mol)
    /// @warning An error is thrown if `material` has elements not present in the chemical system of this material.
    auto addMaterialAmount(Material const& material, double amount) -> void;

    /// Add a given amount of another material in this material.
    /// @param material The other material composing this material.
    /// @param amount The amount of the substance
    /// @param unit The amount unit
    /// @warning An error is thrown if `material` has elements not present in the chemical system of this material.
    auto addMaterialAmount(Material const& material, double amount, Chars unit) -> void;

    /// Add a given mass of another material in this material.
    /// @param material The other material composing this material.
    /// @param mass The mass of the substance
    /// @param unit The mass unit
    /// @warning An error is thrown if `material` has elements not present in the chemical system of this material.
    auto addMaterialMass(Material const& material, double mass, Chars unit) -> void;

    /// Add a given amount or mass of a substance or existing species to the material.
    /// @param substance The chemical formula of the substance or the name of an existing species in the chemical system.
    /// @param value The amount or mass value of the added substance.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if `substance` has elements not present in the chemical system of this material.
    auto add(String const& substance, double value, Chars unit) -> void;

    /// Add a given amount or mass of another material in this material.
    /// @param material The other material composing this material.
    /// @param value The amount or mass value of the added material.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if `material` has elements not present in the chemical system of this material.
    auto add(Material const& material, double value, Chars unit) -> void;

    /// Scale the amount of the material to a desired value.
    /// This method will adjust the amounts of the material's
    /// substances to obtain the desired amount value.
    /// @param value The desired amount of the material
    /// @param unit The amount unit
    auto scaleAmount(double value, Chars unit) -> void;

    /// Scale the mass of the material to a desired value.
    /// This method will adjust the amounts of the material's
    /// substances to obtain the desired mass value.
    /// @param value The desired mass of the material
    /// @param unit The mass unit
    auto scaleMass(double value, Chars unit) -> void;

    /// Scale the amount or mass of the material to a desired value.
    /// This method will adjust the amounts of the material's
    /// substances to obtain the desired amount or mass value.
    /// @param value The desired amount or mass of the material
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    auto scale(double value, Chars unit) -> void;

    /// Return a copy of this material with scaled amount or mass.
    /// @param value The desired amount or mass of the copied material
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    auto with(double value, Chars unit) const -> Material;

    /// Return the underlying chemical system in this material.
    auto system() const -> const ChemicalSystem&;

    /// Return the substances and their amounts in this material.
    auto substances() const -> const Pairs<ChemicalFormula, double>&;

    /// Return the chemical species (as indices) and their amounts in this material.
    auto species() const -> const Pairs<Index, double>&;

    /// Return the accumulated amounts of elements and electric charge in the material.
    auto componentAmounts() const -> ArrayXd;

    /// Return the accumulated amounts of elements in the material.
    auto elementAmounts() const -> ArrayXd;

    /// Return the accumulated electric charge in the material.
    auto charge() const -> double;

    /// Return the amount of the material as the sum of its substance amounts (in mol).
    auto amount() const -> double;

    /// Return the mass of the material as the sum of its substance masses (in kg).
    auto mass() const -> double;

    /// Return the molar mass of the material as the ratio of its mass and amount (in kg/mol).
    auto molarMass() const -> double;

    /// Perform a chemical equilibrium calculation on this material at 25 °C and 1 bar.
    auto equilibrate() -> ChemicalState;

    /// Perform a chemical equilibrium calculation on this material at 25 °C and 1 bar.
    auto equilibrate(const EquilibriumOptions& options) -> ChemicalState;

    /// Perform a chemical equilibrium calculation on this material at 25 °C and 1 bar.
    auto equilibrate(const EquilibriumRestrictions& restrictions) -> ChemicalState;

    /// Perform a chemical equilibrium calculation on this material at 25 °C and 1 bar.
    auto equilibrate(const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options) -> ChemicalState;

    /// Perform a chemical equilibrium calculation on this material at given temperature and pressure.
    auto equilibrate(double T, Chars unitT, double P, Chars unitP) -> ChemicalState;

    /// Perform a chemical equilibrium calculation on this material at given temperature and pressure.
    auto equilibrate(double T, Chars unitT, double P, Chars unitP, const EquilibriumOptions& options) -> ChemicalState;

    /// Perform a chemical equilibrium calculation on this material at given temperature and pressure.
    auto equilibrate(double T, Chars unitT, double P, Chars unitP, const EquilibriumRestrictions& restrictions) -> ChemicalState;

    /// Perform a chemical equilibrium calculation on this material at given temperature and pressure.
    auto equilibrate(double T, Chars unitT, double P, Chars unitP, const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options) -> ChemicalState;

    /// Return the result of the chemical equilibrium calculation performed by `equilibrate`.
    auto result() const -> const EquilibriumResult&;

    /// Return the initial chemical state that is used in the chemical
    /// equilibrium calculation performed by `equilibrate`. By accessing this
    /// initial chemical state, you can see how your material definition was
    /// used to construct an initial guess for the equilibrium calculation.
    /// @param T The temperature for the initial chemical state (in K)
    /// @param P The pressure for the initial chemical state (in Pa)
    auto initialState(double T, double P) const -> ChemicalState;

    /// Return a copy of this material with scaled amount or mass (equivalent to @ref with).
    /// @param value The desired amount or mass of the copied material
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    auto operator()(double value, Chars unit) const -> Material;

private:
    /// The chemical system associated to this material.
    ChemicalSystem m_system;

    /// The substances and their amounts (given as formulas) composing this material.
    Pairs<ChemicalFormula, double> m_substances;

    /// The substances and their amounts (given as existing species) composing this material.
    Pairs<Index, double> m_species;

    /// The result of the chemical equilibrium calculation performed by `equilibrate`
    EquilibriumResult m_result;
};

/// Return a material that is the result of the combination of two others.
auto operator+(const Material& l, const Material& r) -> Material;

/// Output a Material object to an output stream.
auto operator<<(std::ostream& out, const Material& material) -> std::ostream&;

} // namespace Reaktoro
