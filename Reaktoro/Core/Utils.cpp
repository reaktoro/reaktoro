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

#include "Utils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {
namespace detail {

auto molarMasses(const SpeciesList& species) -> ArrayXd
{
    ArrayXd molar_masses(species.size());
    transform(species, molar_masses, [](auto&& s) { return s.molarMass(); });
    return molar_masses;
}

auto computeSpeciesAmount(const ChemicalSystem& system, Index ispecies, real value, const String& unit) -> real
{
    if(unit == "mol")
        return value;

    if(units::convertible(unit, "kg"))
    {
        const auto molarmass = system.species(ispecies).molarMass(); // in kg/mol
        value = units::convert(value, unit, "kg"); // from some mass unit to kg
        return value / molarmass; // from kg to mol
    }
    else return units::convert(value, unit, "mol"); // from some amount unit to mol
}

auto resolveElementIndexAux(const ChemicalSystem& system, Index index) -> Index
{
    return index;
}

auto resolveElementIndexAux(const ChemicalSystem& system, int index) -> Index
{
    return index;
}

auto resolveElementIndexAux(const ChemicalSystem& system, const String& symbol) -> Index
{
    return system.elements().index(symbol);
}

auto resolveElementIndex(const ChemicalSystem& system, StringOrIndex element) -> Index
{
    return std::visit([&](auto&& arg) { return resolveElementIndexAux(system, arg); }, element);
}

auto resolveSpeciesIndexAux(const ChemicalSystem& system, Index index) -> Index
{
    return index;
}

auto resolveSpeciesIndexAux(const ChemicalSystem& system, int index) -> Index
{
    return index;
}

auto resolveSpeciesIndexAux(const ChemicalSystem& system, const String& name) -> Index
{
    return system.species().index(name);
}

auto resolveSpeciesIndex(const ChemicalSystem& system, StringOrIndex species) -> Index
{
return std::visit([&](auto&& arg) { return resolveSpeciesIndexAux(system, arg); }, species);
}

auto resolvePhaseIndexAux(const ChemicalSystem& system, Index index) -> Index
{
    return index;
}

auto resolvePhaseIndexAux(const ChemicalSystem& system, int index) -> Index
{
    return index;
}

auto resolvePhaseIndexAux(const ChemicalSystem& system, const String& name) -> Index
{
    return system.phases().index(name);
}

auto resolvePhaseIndex(const ChemicalSystem& system, StringOrIndex phase) -> Index
{
    return std::visit([&](auto&& arg) { return resolvePhaseIndexAux(system, arg); }, phase);
}

} // namespace detail
} // namespace Reaktoro
