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
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

/// The class used to define reactivity restrictions in a chemical equilibrium calculation.
class EquilibriumRestrictions
{
public:
    /// Construct an EquilibriumRestrictions object.
    explicit EquilibriumRestrictions(const ChemicalSystem& system);

    /// Prevent the amount of a species from changing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species whose amount cannot change.
    auto cannotReact(Index ispecies) -> void;

    /// Prevent the amount of a species from changing during the chemical equilibrium calculation.
    /// @param species The name of the species whose amount cannot change.
    auto cannotReact(String species) -> void;

    /// Prevent the amount of a species from increasing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species.
    auto cannotIncrease(Index ispecies) -> void;

    /// Prevent the amount of a species from increasing during the chemical equilibrium calculation.
    /// @param species The name of the species.
    auto cannotIncrease(String species) -> void;

    /// Prevent the amount/mass of a species from increasing above a value during the chemical equilibrium calculation.
    /// @param ispecies The index of the species.
    /// @param value The upper bound value of the species amount/mass.
    /// @param unit The unit of the upper bound value (must be convertible to mol or kg).
    auto cannotIncreaseAbove(Index ispecies, double value, String unit="mol") -> void;

    /// Prevent the amount/mass of a species from increasing above a value during the chemical equilibrium calculation.
    /// @param species The name of the species.
    /// @param value The upper bound value of the species amount/mass.
    /// @param unit The unit of the upper bound value (must be convertible to mol or kg).
    auto cannotIncreaseAbove(String species, double value, String unit="mol") -> void;

    /// Prevent the amount of a species from decreasing during the chemical equilibrium calculation.
    /// @param ispecies The index of the species.
    auto cannotDecrease(Index ispecies) -> void;

    /// Prevent the amount of a species from decreasing during the chemical equilibrium calculation.
    /// @param species The name of the species.
    auto cannotDecrease(String species) -> void;

    /// Prevent the amount/mass of a species from decreasing below a value during the chemical equilibrium calculation.
    /// @param ispecies The index of the species.
    /// @param value The lower bound value of the species amount/mass.
    /// @param unit The unit of the lower bound value (must be convertible to mol or kg).
    auto cannotDecreaseBelow(Index ispecies, double value, String unit="mol") -> void;

    /// Prevent the amount/mass of a species from decreasing below a value during the chemical equilibrium calculation.
    /// @param species The name of the species.
    /// @param value The lower bound value of the species amount/mass.
    /// @param unit The unit of the lower bound value (must be convertible to mol or kg).
    auto cannotDecreaseBelow(String species, double value, String unit="mol") -> void;

    /// Allow the amount of a species to change during the chemical equilibrium calculation.
    /// @param ispecies The index of the species whose amount is allowed to change.
    auto canReact(Index ispecies) -> void;

    /// Allow the amount of a species to change during the chemical equilibrium calculation.
    /// @param species The name of the species whose amount is allowed to change.
    auto canReact(String species) -> void;

    /// Allow the amount of a species to increase during the chemical equilibrium calculation.
    /// @param ispecies The index of the species whose amount is allowed to increase.
    auto canIncrease(Index ispecies) -> void;

    /// Allow the amount of a species to increase during the chemical equilibrium calculation.
    /// @param species The name of the species whose amount is allowed to increase.
    auto canIncrease(String species) -> void;

    /// Allow the amount of a species to decrease during the chemical equilibrium calculation.
    /// @param ispecies The index of the species whose amount is allowed to decrease.
    auto canDecrease(Index ispecies) -> void;

    /// Allow the amount of a species to decrease during the chemical equilibrium calculation.
    /// @param species The name of the species whose amount is allowed to decrease.
    auto canDecrease(String species) -> void;

    /// Return the chemical system associated with the equilibrium conditions.
    auto system() const -> const ChemicalSystem&;

    /// Return the indices of the species whose amounts cannot increase.
    auto indicesSpeciesCannotIncrease() const -> Set<Index> const&;

    /// Return the indices of the species whose amounts cannot decrease.
    auto indicesSpeciesCannotDecrease() const -> Set<Index> const&;

    /// Return the indices of the species whose amounts cannot increase above a given value.
    auto indicesSpeciesCannotIncreaseAbove() const -> Map<Index, double> const&;

    /// Return the indices of the species whose amounts cannot decrease below a given value.
    auto indicesSpeciesCannotDecreaseBelow() const -> Map<Index, double> const&;

private:
    /// The chemical system associated with the equilibrium restrictions.
    ChemicalSystem msystem;

    /// The indices of the species whose amounts cannot increase above its current value.
    Set<Index> species_cannot_increase;

    /// The indices of the species whose amounts cannot decrease below its current value.
    Set<Index> species_cannot_decrease;

    /// The indices and upper bound amounts (in mol) of the species with given upper bounds.
    Map<Index, double> species_cannot_increase_above;

    /// The indices and lower bound amounts (in mol) of the species with given lower bounds.
    Map<Index, double> species_cannot_decrease_below;
};

} // namespace Reaktoro
