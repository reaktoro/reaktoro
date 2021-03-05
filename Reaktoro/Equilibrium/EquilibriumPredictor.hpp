// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class EquilibriumConditions;
class EquilibriumSensitivity;

/// A class that still needs to be commented.
class EquilibriumPredictor
{
public:
    /// Construct a EquilibriumPredictor object.
    /// @param state0 The reference chemical equilibrium state from which first-order Taylor predictions are made.
    /// @param sensitivity0 The sensitivity derivatives of the chemical equilibrium state at the reference point.
    EquilibriumPredictor(const ChemicalState& state0, const EquilibriumSensitivity& sensitivity0);

    /// Construct a copy of a EquilibriumPredictor object.
    EquilibriumPredictor(const EquilibriumPredictor& other);

    /// Destroy this EquilibriumPredictor object.
    ~EquilibriumPredictor();

    /// Assign a copy of a EquilibriumPredictor object to this.
    auto operator=(EquilibriumPredictor other) -> EquilibriumPredictor&;

    /// Perform a first-order Taylor prediction of the chemical state at given conditions.
    auto predict(ChemicalState& state, const EquilibriumConditions& conditions) -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
