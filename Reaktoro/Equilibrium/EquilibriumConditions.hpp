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
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>

namespace Reaktoro {

/// The class used to define conditions to be satisfied at chemical equilibrium.
class EquilibriumConditions
{
public:
    /// Construct an EquilibriumConditions object.
    explicit EquilibriumConditions(const EquilibriumSpecs& specs);

    //=================================================================================================
    //
    // METHODS TO SPECIFY THERMODYNAMIC CONDITIONS
    //
    //=================================================================================================

    /// Specify the **temperature** of the system at chemical equilibrium.
    /// @param value The temperature of the system
    /// @param unit The unit of the temperature value (must be convertible to K)
    auto temperature(real value, String unit="K") -> void;

    /// Specify the **pressure** of the system at chemical equilibrium.
    /// @param value The pressure of the system
    /// @param unit The unit of the pressure value (must be convertible to Pa)
    auto pressure(real value, String unit="Pa") -> void;

    /// Specify the **volume** of the system at chemical equilibrium.
    /// @param value The volume of the system
    /// @param unit The unit of the volume value (must be convertible to m@sup{3})
    auto volume(real value, String unit="m3") -> void;

    /// Specify the **internal energy** of the system at chemical equilibrium.
    /// @param value The internal energy of the system
    /// @param unit The unit of the internal energy value (must be convertible to J)
    auto internalEnergy(real value, String unit="J") -> void;

    /// Specify the **enthalpy** of the system at chemical equilibrium.
    /// @param value The enthalpy of the system
    /// @param unit The unit of the enthalpy value (must be convertible to J)
    auto enthalpy(real value, String unit="J") -> void;

    /// Specify the **Gibbs energy** of the system at chemical equilibrium.
    /// @param value The Gibbs energy of the system
    /// @param unit The unit of the Gibbs energy value (must be convertible to J)
    auto gibbsEnergy(real value, String unit="J") -> void;

    /// Specify the **Helmholtz energy** of the system at chemical equilibrium.
    /// @param value The Helmholtz energy of the system
    /// @param unit The unit of the Helmholtz energy value (must be convertible to J)
    auto helmholtzEnergy(real value, String unit="J") -> void;

    /// Specify the **entropy** of the system at chemical equilibrium.
    /// @param value The entropy of the system
    /// @param unit The unit of the entropy value (must be convertible to J/K)
    auto entropy(real value, String unit="J/K") -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY CHEMICAL POTENTIAL CONSTRAINTS
    //
    //=================================================================================================

    /// Specify the **chemical potential** of a substance at chemical equilibrium.
    /// @param substance The chemical formula of the substance.
    /// @param value The constrained chemical potential value.
    /// @param unit The unit for the constrained chemical potential value (must be convertible to J/mol).
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a chemical potential constraint for the substance.
    auto chemicalPotential(String substance, real value, String unit="J/mol") -> void;

    /// Specify the **ln activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species.
    /// @param value The constrained ln activity value.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider an activity constraint for the species.
    auto lnActivity(String species, real value) -> void;

    /// Specify the **lg activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species.
    /// @param value The constrained lg activity value.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider an activity constraint for the species.
    auto lgActivity(String species, real value) -> void;

    /// Specify the **activity** of a species at chemical equilibrium.
    /// @param species The name of the chemical species.
    /// @param value The constrained activity value.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider an activity constraint for the species.
    auto activity(String species, real value) -> void;

    /// Specify the **fugacity** of a gaseous species at chemical equilibrium.
    /// @param species The name of the gaseous species.
    /// @param value The constrained fugacity value.
    /// @param unit The unit for the constrained fugacity value (must be convertible to Pa).
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a fugacity constraint for the gas.
    auto fugacity(String species, real value, String unit) -> void;

    /// Specify the *pH* at chemical equilibrium.
    /// @param value The constrained value for pH.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a pH constraint.
    auto pH(real value) -> void;

    /// Specify the *pMg* at chemical equilibrium.
    /// @param value The constrained value for pMg.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a pMg constraint.
    auto pMg(real value) -> void;

    /// Specify the *pE* at chemical equilibrium.
    /// @param value The constrained value for pE.
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider a pE constraint.
    auto pE(real value) -> void;

    /// Specify the *Eh* at chemical equilibrium.
    /// @param value The constrained value for Eh.
    /// @param unit The unit of the constrained value for Eh (must be convertible to V).
    /// @warning An error is thrown if the specifications for the chemical equilibrium calculation
    /// @warning do not consider an Eh constraint.
    auto Eh(real value, String unit="V") -> void;

    //=================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================

    /// Set the initial amounts of the conservative components in the equilibrium calculation.
    /// These amounts are conserved at chemical equilibrium if the system is
    /// closed. If the system is open to one or more substances, these initial
    /// amounts will differ from those computed at chemical equilibrium. The
    /// difference correspond to how much each titrant (i.e., the substance for
    /// which the system is open to) entered/leaved the system.
    /// @param values The initial amounts of the components in the system.
    auto initialComponentAmounts(VectorXrConstRef values) -> void;

    /// Return the initial amounts of the conservative components in the equilibrium calculation.
    auto initialComponentAmounts() const -> VectorXrConstRef;

    /// Return the chemical system associated with the equilibrium conditions.
    auto system() const -> const ChemicalSystem&;

    /// Return the input parameters associated with the equilibrium conditions.
    auto params() const -> Params;

private:
    /// The chemical system associated with the given equilibrium specifications.
    const ChemicalSystem m_system;

    /// The input parameters defining the equilibrium conditions according to given equilibrium specifications.
    Params m_params;

    /// The initial amounts of the conservative components in the equilibrium calculation.
    VectorXr m_initial_component_amounts;
};

} // namespace Reaktoro
