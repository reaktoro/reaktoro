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
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Common/Algorithms.hpp>

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

auto resolveElementIndexAux(const Phase& phase, Index index) -> Index
{
    return index;
}

auto resolveElementIndexAux(const Phase& phase, int index) -> Index
{
    return index;
}

auto resolveElementIndexAux(const Phase& phase, const String& symbol) -> Index
{
    return phase.elements().index(symbol);
}

auto resolveElementIndex(const Phase& phase, StringOrIndex element) -> Index
{
    return std::visit([&](auto&& arg) { return resolveElementIndexAux(phase, arg); }, element);
}

auto resolveSpeciesIndexAux(const Phase& phase, Index index) -> Index
{
    return index;
}

auto resolveSpeciesIndexAux(const Phase& phase, int index) -> Index
{
    return index;
}

auto resolveSpeciesIndexAux(const Phase& phase, const String& name) -> Index
{
    return phase.species().index(name);
}

auto resolveSpeciesIndex(const Phase& phase, StringOrIndex element) -> Index
{
    return std::visit([&](auto&& arg) { return resolveSpeciesIndexAux(phase, arg); }, element);
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

auto assembleFormulaMatrix(const SpeciesList& species, const ElementList& elements) -> MatrixXd
{
    const auto num_elements = elements.size();
    const auto num_components = num_elements + 1;
    const auto num_species = species.size();
    MatrixXd A(num_components, num_species);
    for(auto i = 0; i < num_species; ++i)
        for(auto j = 0; j < num_elements; ++j)
            A(j, i) = species[i].elements().coefficient(elements[j].symbol());
    for(auto i = 0; i < num_species; ++i)
        A(num_elements, i) = species[i].charge();
    return A;
}

auto extractNames(const SpeciesList& list) -> Strings
{
    return vectorize(list, RKT_LAMBDA(x, x.name()));
}

auto extractNames(const ElementList& list) -> Strings
{
    return vectorize(list, RKT_LAMBDA(x, x.name()));
}

auto extractNames(const PhaseList& list) -> Strings
{
    return vectorize(list, RKT_LAMBDA(x, x.name()));
}

} // namespace detail
} // namespace Reaktoro
