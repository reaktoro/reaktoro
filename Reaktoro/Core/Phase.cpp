// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "Phase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/PhaseChemicalProperties.hpp>
#include <Reaktoro/Core/PhaseThermoProperties.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

struct Phase::Impl
{
    /// The name of the phase
    std::string name;

    /// The type of the phase.
    PhaseType type = PhaseType::Solid;

    /// The list of Species instances defining the phase
    std::vector<Species> species;

    /// The list of Element instances in the phase
    std::vector<Element> elements;

    /// The function that calculates the standard thermodynamic properties of the phase and its species
    PhaseThermoModel thermo_model;

    /// The function that calculates the chemical properties of the phase and its species
    PhaseChemicalModel chemical_model;

    // The molar masses of the species
    Vector molar_masses;
};

Phase::Phase()
: pimpl(new Impl())
{}

auto Phase::setName(std::string name) -> void
{
    pimpl->name = name;
}

auto Phase::setType(PhaseType type) -> void
{
    pimpl->type = type;
}

auto Phase::setSpecies(const std::vector<Species>& species) -> void
{
    pimpl->species = species;
    pimpl->molar_masses = molarMasses(species);
}

auto Phase::setThermoModel(const PhaseThermoModel& model) -> void
{
    pimpl->thermo_model = model;
}

auto Phase::setChemicalModel(const PhaseChemicalModel& model) -> void
{
    pimpl->chemical_model = model;
}

auto Phase::numElements() const -> unsigned
{
    return elements().size();
}

auto Phase::numSpecies() const -> unsigned
{
    return species().size();
}

auto Phase::name() const -> std::string
{
    return pimpl->name;
}

auto Phase::type() const -> PhaseType
{
    return pimpl->type;
}

auto Phase::elements() const -> const std::vector<Element>&
{
    return pimpl->elements;
}

auto Phase::elements() -> std::vector<Element>&
{
    return pimpl->elements;
}

auto Phase::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto Phase::species() -> std::vector<Species>&
{
    return pimpl->species;
}

auto Phase::species(Index index) const -> const Species&
{
    return pimpl->species[index];
}

auto Phase::isFluid() const -> bool
{
    return !isSolid();
}

auto Phase::isSolid() const -> bool
{
    return type() == PhaseType::Solid;
}

auto Phase::thermoModel() const -> const PhaseThermoModel&
{
    return pimpl->thermo_model;
}

auto Phase::chemicalModel() const -> const PhaseChemicalModel&
{
    return pimpl->chemical_model;
}

auto Phase::indexSpecies(std::string name) const -> Index
{
    return index(name, species());
}

auto Phase::indexSpeciesWithError(std::string name) const -> Index
{
    const Index index = indexSpecies(name);
    Assert(index < numSpecies(),
        "Could not get the index of species `" + name + "`.",
        "There is no species called `" + name + "` in the phase.");
    return index;
}

auto Phase::indexSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    return indexAny(names, species());
}

auto Phase::indexSpeciesAnyWithError(const std::vector<std::string>& names) const -> Index
{
    const Index index = indexSpeciesAny(names);
    Assert(index < numSpecies(),
        "Could not get the index of the species with "
        "any of the following names `" + join(names, ", ") + "`.",
        "There is no species in phase `" + name() + "` with any of these names.");
    return index;
}

auto Phase::properties(double T, double P) const -> PhaseThermoProperties
{
    PhaseThermoProperties prop(*this);
    prop.update(T, P);
    return prop;
}

auto Phase::properties(double T, double P, const Vector& n) const -> PhaseChemicalProperties
{
    PhaseChemicalProperties prop(*this);
    prop.update(T, P, n);
    return prop;
}

auto operator<(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
