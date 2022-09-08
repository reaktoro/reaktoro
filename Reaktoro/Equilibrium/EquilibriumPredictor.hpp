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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class EquilibriumConditions;
class EquilibriumSensitivity;

/// Used to predict a chemical equilibrium state at given conditions using first-order Taylor approximation.
class EquilibriumPredictor
{
public:
    /// Construct a EquilibriumPredictor object.
    /// @param state0 The reference chemical equilibrium state from which first-order Taylor predictions are made.
    /// @param sensitivity0 The sensitivity derivatives of the chemical equilibrium state at the reference point.
    EquilibriumPredictor(ChemicalState const& state0, EquilibriumSensitivity const& sensitivity0);

    /// Construct a copy of a EquilibriumPredictor object.
    EquilibriumPredictor(EquilibriumPredictor const& other);

    /// Destroy this EquilibriumPredictor object.
    ~EquilibriumPredictor();

    /// Assign a copy of a EquilibriumPredictor object to this.
    auto operator=(EquilibriumPredictor other) -> EquilibriumPredictor&;

    /// Perform a first-order Taylor prediction of the chemical state at given conditions.
    /// @param[out] state The predicted chemical equilibrium state
    /// @param conditions The conditons at which the chemical equilibrium state must be satisfied
    auto predict(ChemicalState& state, EquilibriumConditions const& conditions) const -> void;

    /// Perform a first-order Taylor prediction of the chemical state at given conditions.
    /// @param[out] state The predicted chemical equilibrium state
    /// @param dw The change in the values of the input variables *w*.
    /// @param dc The change in the values of the initial amounts of conservative components *c*.
    auto predict(ChemicalState& state, VectorXdConstRef const& dw, VectorXdConstRef const& dc) const -> void;

    /// Perform a first-order Taylor prediction of the chemical potential of a species at given conditions.
    auto speciesChemicalPotentialPredicted(Index ispecies, VectorXdConstRef const& dw, VectorXdConstRef const& dc) const -> double;

    /// Return the chemical potential of a species at given reference conditions.
    auto speciesChemicalPotentialReference(Index ispecies) const -> double;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
