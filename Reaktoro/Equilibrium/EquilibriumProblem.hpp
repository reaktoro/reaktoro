// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;

/// The class used to define chemical equilibrium problems.
class EquilibriumProblem : public EquilibriumConditions, public EquilibriumRestrictions
{
public:
    /// Construct an EquilibriumProblem object with given chemical system.
    explicit EquilibriumProblem(const ChemicalSystem& system);

    /// Construct an EquilibriumProblem object with given specifications.
    explicit EquilibriumProblem(const EquilibriumSpecs& specs);

    /// Specify an initial temperature condition for the chemical system (in K).
    /// @param value The initial temperature of the system.
    auto startWithTemperature(real value) -> void;

    /// Specify an initial temperature condition for the chemical system.
    /// @param value The initial temperature of the system.
    /// @param unit The temperature unit (must be convertible to K).
    auto startWithTemperature(real value, String unit) -> void;

    /// Specify an initial pressure condition for the chemical system.
    /// @param value The initial pressure of the system.
    auto startWithPressure(real value) -> void;

    /// Specify an initial pressure condition for the chemical system.
    /// @param value The initial pressure of the system.
    /// @param unit The pressure unit (must be convertible to Pa).
    auto startWithPressure(real value, String unit) -> void;

    /// Specify an initial condition for the species amounts.
    /// @param n The array with the initial amounts of the species (in mol).
    auto startWithSpeciesAmounts(ArrayXrConstRef n) -> void;

    /// @param species The name of the species in the chemical system.
    /// @param value The amount or mass value of the species.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    /// @warning An error is thrown if the chemical system has no species with name @p species.
    auto startWith(String species, real value, String unit) -> void;

    /// Specify an initial condition for the abundance of a chemical species.
    /// @param ispecies The index of the species in the chemical system.
    /// @param value The amount or mass value of the species.
    /// @param unit The amount or mass unit (must be convertible to mol or kg).
    auto startWith(Index ispecies, real value, String unit) -> void;

    /// Specify an initial condition for temperature, pressure, and species amounts with a given chemical state.
    /// @note This method overwrites all previous `startWith`, `startWithTemperature`,
    /// @note `startWithPressure`, and `startWithSpeciesAmounts` method calls.
    /// @param state The chemical state representing the initial conditions of the system.
    auto startWithState(const ChemicalState& state) -> void;

    /// Specify the initial condition for the amounts of the conservative components.
    /// These component amounts are conserved at chemical equilibrium only if
    /// the system is closed. If the system is open to one or more substances,
    /// these given initial component amounts will differ from those computed
    /// at chemical equilibrium. The difference correspond to how much each
    /// titrant (i.e., the substance for which the system is open to) entered
    /// or leaved the system.
    ///
    /// @note This method has higher priority than the other methods that
    /// specify initial conditions for the abundance of chemical species:
    ///
    ///    * @ref EquilibriumProblem::startWith
    ///    * @ref EquilibriumProblem::startWithSpeciesAmounts
    ///    * @ref EquilibriumProblem::startWithState
    ///
    /// By using this method, the chemical equilibrium calculation will use the
    /// specified component amounts instead of computing them from the initial
    /// condition for the species amounts. In this case, the provided
    /// conditions for species amounts will be solely used as initial guess for
    /// the equilibrium computation. The provided species amounts may not be
    /// consistent with the given component amounts here (i.e., there is mass
    /// imbalance between initial condition for species amounts and the amounts
    /// of the components that constitute the species). At the end of a
    /// successfull calculation, this imbalance is resolved, and the species
    /// amounts are consistent with the component amounts in terms of mass
    /// conservation. This is true provided the chemical system is not open to
    /// other substances, in which case, mass conservation also considers
    /// these substances as well.
    ///
    /// @param b The array with the initial amounts of the components (in mol)
    auto startWithComponentAmounts(ArrayXrConstRef b) -> void;

    /// Return the initial temperature of the system in the equilibrium calculation (in K).
    auto initialTemperature() const -> real;

    /// Return the initial pressure of the system in the equilibrium calculation (in Pa).
    auto initialPressure() const -> real;

    /// Return the initial amounts of the species in the equilibrium calculation (in mol).
    auto initialSpeciesAmounts() const -> ArrayXrConstRef;

    /// Return the initial amounts of the conservative components in the equilibrium calculation (in mol).
    auto initialComponentAmounts() const -> ArrayXr;

    using EquilibriumConditions::system;

private:

    /// The initial temperature of the system in the equilibrium calculation (in K).
    real m_initial_temperature = NaN;

    /// The initial pressure of the system in the equilibrium calculation (in Pa).
    real m_initial_pressure = NaN;

    /// The initial amounts of the species in the equilibrium calculation (in mol).
    ArrayXr m_initial_species_amounts;

    /// The initial amounts of the conservative components in the equilibrium calculation (in mol).
    ArrayXr m_initial_component_amounts;

};

} // namespace Reaktoro
