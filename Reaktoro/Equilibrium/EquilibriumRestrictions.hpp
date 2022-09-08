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
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

/// The class used to define reactivity restrictions in a chemical equilibrium calculation.
class EquilibriumRestrictions
{
public:
    /// Construct an EquilibriumRestrictions object.
    explicit EquilibriumRestrictions(ChemicalSystem const& system);

    //=================================================================================================
    //
    // METHODS THAT INTRODUCE RESTRICTIONS ON HOW MUCH SPECIES CAN REACT
    //
    //=================================================================================================

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
    auto cannotIncreaseAbove(Index ispecies, double value, Chars unit="mol") -> void;

    /// Prevent the amount/mass of a species from increasing above a value during the chemical equilibrium calculation.
    /// @param species The name of the species.
    /// @param value The upper bound value of the species amount/mass.
    /// @param unit The unit of the upper bound value (must be convertible to mol or kg).
    auto cannotIncreaseAbove(String species, double value, Chars unit="mol") -> void;

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
    auto cannotDecreaseBelow(Index ispecies, double value, Chars unit="mol") -> void;

    /// Prevent the amount/mass of a species from decreasing below a value during the chemical equilibrium calculation.
    /// @param species The name of the species.
    /// @param value The lower bound value of the species amount/mass.
    /// @param unit The unit of the lower bound value (must be convertible to mol or kg).
    auto cannotDecreaseBelow(String species, double value, Chars unit="mol") -> void;

    //=================================================================================================
    //
    // METHODS THAT REMOVE EXISTING REACTIVITY RESTRICTIONS ON SPECIES
    //
    //=================================================================================================

    /// Remove any previously set increase/decrease restriction on the amount of a species.
    /// @param ispecies The index of the species whose amount is now allowed to change without bounds.
    auto canReactFreely(Index ispecies) -> void;

    /// Remove any previously set increase/decrease restriction on the amount of a species.
    /// @param species The name of the species whose amount is now allowed to change without bounds.
    auto canReactFreely(String species) -> void;

    /// Remove any previously set increase restriction on the amount of a species.
    /// @param ispecies The index of the species whose amount is now allowed to increase without bounds.
    auto canIncreaseFreely(Index ispecies) -> void;

    /// Remove any previously set increase restriction on the amount of a species.
    /// @param species The name of the species whose amount is now allowed to increase without bounds.
    auto canIncreaseFreely(String species) -> void;

    /// Remove any previously set decrease restriction on the amount of a species.
    /// @param ispecies The index of the species whose amount is now allowed to decrease without bounds.
    auto canDecreaseFreely(Index ispecies) -> void;

    /// Remove any previously set decrease restriction on the amount of a species.
    /// @param species The name of the species whose amount is now allowed to decrease without bounds.
    auto canDecreaseFreely(String species) -> void;

    //=================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================

    /// Return the chemical system associated with the equilibrium conditions.
    auto system() const -> ChemicalSystem const&;

    /// Return the indices of the species whose amounts cannot increase.
    auto speciesCannotIncrease() const -> Set<Index> const&;

    /// Return the indices of the species whose amounts cannot decrease.
    auto speciesCannotDecrease() const -> Set<Index> const&;

    /// Return the indices of the species whose amounts cannot increase above a given value.
    auto speciesCannotIncreaseAbove() const -> Map<Index, double> const&;

    /// Return the indices of the species whose amounts cannot decrease below a given value.
    auto speciesCannotDecreaseBelow() const -> Map<Index, double> const&;

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
