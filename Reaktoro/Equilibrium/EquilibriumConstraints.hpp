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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class ChemicalSystem;

/// The class used to define equilibrium constraints.
class EquilibriumConstraints
{
public:
    /// Construct a default EquilibriumConstraints object.
    EquilibriumConstraints();

    /// Enforce a value for the **temperature** of the system at the chemical equilibrium state.
    /// @param value The constrained temperature of the system
    /// @param unit The unit of the constrained temperature value (must be convertible to K)
    auto temperature(real value, String unit) -> void;

    /// Enforce a value for the **pressure** of the system at the chemical equilibrium state.
    /// @param value The constrained pressure of the system
    /// @param unit The unit of the constrained pressure value (must be convertible to Pa)
    auto pressure(real value, String unit) -> void;

    /// Enforce a value for the **volume** of the system at the chemical equilibrium state.
    /// @param value The constrained volume of the system
    /// @param unit The unit of the constrained volume value (must be convertible to m@sup{3})
    auto volume(real value, String unit) -> void;

    /// Enforce a value for the **internal energy** of the system at the chemical equilibrium state.
    /// @param value The constrained internal energy of the system
    /// @param unit The unit of the constrained internal energy value (must be convertible to J)
    auto internalEnergy(real value, String unit) -> void;

    /// Enforce a value for the **enthalpy** of the system at the chemical equilibrium state.
    /// @param value The constrained enthalpy of the system
    /// @param unit The unit of the constrained enthalpy value (must be convertible to J)
    auto enthalpy(real value, String unit) -> void;

    /// Enforce a value for the **Gibbs energy** of the system at the chemical equilibrium state.
    /// @param value The constrained Gibbs energy of the system
    /// @param unit The unit of the constrained Gibbs energy value (must be convertible to J)
    auto gibbsEnergy(real value, String unit) -> void;

    /// Enforce a value for the **Helmholtz energy** of the system at the chemical equilibrium state.
    /// @param value The constrained Helmholtz energy of the system
    /// @param unit The unit of the constrained Helmholtz energy value (must be convertible to J)
    auto helmholtzEnergy(real value, String unit) -> void;

    /// Enforce a value for the **entropy** of the system at the chemical equilibrium state.
    /// @param value The constrained entropy of the system
    /// @param unit The unit of the constrained entropy value (must be convertible to J/K)
    auto entropy(real value, String unit) -> void;

    /// Enforce a value for the **chemical potential** of a species at the chemical equilibrium state.
    /// @param species The name of the chemical species as found in the database under use
    /// @param value The constrained chemical potential value
    /// @param unit The unit for the constrained chemical potential value (must be convertible to J/mol)
    /// @note The chemical species need not be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto chemicalPotential(String species, real value, String unit) -> void;

    /// Enforce a value for the **chemical potential** of a species at the chemical equilibrium state.
    /// @param species The name of the chemical species as found in the database under use
    /// @param fn The function that computes the constrained chemical potential value at given temperature and pressure
    /// @param unit The unit for the constrained chemical potential value (must be convertible to J/mol)
    /// @note The chemical species need not be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto chemicalPotential(String species, Fn<real(real,real)> fn, String unit) -> void;

    /// Enforce a value for the **activity** of a species at the chemical equilibrium state.
    /// @param species The name of the chemical species as found in the database under use
    /// @param value The constrained activity value
    /// @note The chemical species need not be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    auto activity(String species, real value) -> void;

    /// Enforce a value for the **fugacity** of a gaseous species at the chemical equilibrium state.
    /// @param species The name of the gaseous species as found in the database under use
    /// @param value The constrained fugacity value
    /// @param unit The unit for the constrained fugacity value (must be convertible to Pa)
    /// @note The gaseous species need not be in the chemical system; only in the database.
    /// @warning An error will be raised if the database does not contain a species with given name.
    /// @warning It is also not permitted to use a species that is not a gas.
    auto fugacity(String species, real value, String unit) -> void;

    /// Enforce a value for pH at the chemical equilibrium state.
    /// @param value The constrained value for pH
    auto pH(real value) -> void;

    /// Enforce a value for pe at the chemical equilibrium state.
    /// @param value The constrained value for pe
    auto pe(real value) -> void;

    /// Enforce a value for Eh at the chemical equilibrium state.
    /// @param value The constrained value for Eh
    /// @param unit The unit of the constrained value for Eh (must be convertible to V)
    auto Eh(real value, String unit) -> void;
};


} // namespace Reaktoro
