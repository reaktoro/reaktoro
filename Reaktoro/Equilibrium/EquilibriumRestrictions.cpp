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

#include "EquilibriumRestrictions.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>

namespace Reaktoro {
namespace {

/// Compute the amount of a species given a value in mass or in amount.
auto computeSpeciesAmount(const ChemicalSystem& system, Index ispecies, double value, const String& unit)
{
    if(unit == "mol")
        return value;

    if(units::convertible(unit, "kg"))
    {
        const auto molarmass = system.species(ispecies).molarMass(); // in mol/kg
        value = units::convert(value, unit, "kg"); // from some mass unit to kg
        return value * molarmass; // from kg to mol
    }
    else return units::convert(value, unit, "mol"); // from some amount unit to mol
}

} // namespace

EquilibriumRestrictions::EquilibriumRestrictions(const ChemicalSystem& system)
: msystem(system)
{}

auto EquilibriumRestrictions::cannotReact(Index ispecies) -> void
{
    cannotIncrease(ispecies);
    cannotDecrease(ispecies);
}

auto EquilibriumRestrictions::cannotReact(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotReact(ispecies);
}

auto EquilibriumRestrictions::cannotIncrease(Index ispecies) -> void
{
    const auto numspecies = system().species().size();
    error(ispecies < numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    species_cannot_increase.insert(ispecies);
}

auto EquilibriumRestrictions::cannotIncrease(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotIncrease(ispecies);
}

auto EquilibriumRestrictions::cannotIncreaseAbove(Index ispecies, double value, String unit) -> void
{
    const auto numspecies = system().species().size();
    error(ispecies < numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    value = computeSpeciesAmount(msystem, ispecies, value, unit);
    species_cannot_increase_above.insert_or_assign(ispecies, value);
}

auto EquilibriumRestrictions::cannotIncreaseAbove(String species, double value, String unit) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotIncreaseAbove(ispecies, value, unit);
}

auto EquilibriumRestrictions::cannotDecrease(Index ispecies) -> void
{
    const auto numspecies = system().species().size();
    error(ispecies < numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    species_cannot_decrease.insert(ispecies);
}

auto EquilibriumRestrictions::cannotDecrease(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotDecrease(ispecies);
}

auto EquilibriumRestrictions::cannotDecreaseBelow(Index ispecies, double value, String unit) -> void
{
    const auto numspecies = system().species().size();
    error(ispecies < numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    value = computeSpeciesAmount(msystem, ispecies, value, unit);
    species_cannot_decrease_below.insert_or_assign(ispecies, value);
}

auto EquilibriumRestrictions::cannotDecreaseBelow(String species, double value, String unit) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    cannotDecreaseBelow(ispecies, value, unit);
}

auto EquilibriumRestrictions::canReact(Index ispecies) -> void
{
    canIncrease(ispecies);
    canDecrease(ispecies);
}

auto EquilibriumRestrictions::canReact(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    canReact(ispecies);
}

auto EquilibriumRestrictions::canIncrease(Index ispecies) -> void
{
    const auto numspecies = system().species().size();
    error(ispecies < numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    species_cannot_increase.erase(ispecies);
}

auto EquilibriumRestrictions::canIncrease(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    species_cannot_increase.erase(ispecies);
}

auto EquilibriumRestrictions::canDecrease(Index ispecies) -> void
{
    const auto numspecies = system().species().size();
    error(ispecies < numspecies, "Given species index `", ispecies, "` is out of bounds (number of species is ", numspecies, ").");
    species_cannot_decrease.erase(ispecies);
}

auto EquilibriumRestrictions::canDecrease(String species) -> void
{
    const auto ispecies = system().species().indexWithName(species);
    species_cannot_decrease.erase(ispecies);
}

auto EquilibriumRestrictions::system() const -> const ChemicalSystem&
{
    return msystem;
}

auto EquilibriumRestrictions::indicesSpeciesCannotIncrease() const -> Set<Index> const&
{
    return species_cannot_increase;
}

auto EquilibriumRestrictions::indicesSpeciesCannotDecrease() const -> Set<Index> const&
{
    return species_cannot_decrease;
}

auto EquilibriumRestrictions::indicesSpeciesCannotIncreaseAbove() const -> Map<Index, double> const&
{
    return species_cannot_increase_above;
}

auto EquilibriumRestrictions::indicesSpeciesCannotDecreaseBelow() const -> Map<Index, double> const&
{
    return species_cannot_decrease_below;
}

} // namespace Reaktoro
