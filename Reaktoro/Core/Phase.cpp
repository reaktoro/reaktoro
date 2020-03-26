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
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

struct Phase::Impl
{
    /// The name of the phase
    std::string name;

    /// The type of the phase.
    std::string type;

    /// The physical state of the phase.
    PhasePhysicalState state = PhasePhysicalState::Solid;

    /// The list of Species instances defining the phase
    std::vector<Species> species;

    /// The list of Element instances in the phase
    std::vector<Element> elements;

    /// The standard thermodynamic model of the phase.
    PhaseStandardThermoModelFn standard_thermo_model_fn;

    /// The activity model of the phase.
    PhaseActivityModelFn activity_model_fn;
};

Phase::Phase()
: pimpl(new Impl())
{}

Phase::Phase(const Phase& other)
: pimpl(new Impl(*other.pimpl))
{}

Phase::~Phase()
{}

auto Phase::operator=(Phase other) -> Phase&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Phase::setName(std::string name) -> void
{
    pimpl->name = name;
}

auto Phase::setType(std::string type) -> void
{
    pimpl->type = type;
}

auto Phase::setPhysicalState(PhasePhysicalState state) -> void
{
    pimpl->state = state;
}

auto Phase::setSpecies(const std::vector<Species>& species) -> void
{
    pimpl->species = species;
}

auto Phase::setStandardThermoModel(const PhaseStandardThermoModelFn& model) -> void
{
    pimpl->standard_thermo_model_fn = model;
}

auto Phase::setActivityModel(const PhaseActivityModelFn& model) -> void
{
    pimpl->activity_model_fn = model;
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

auto Phase::type() const -> std::string
{
    return pimpl->type;
}

auto Phase::physicalState() const -> PhasePhysicalState
{
    return pimpl->state;
}

auto Phase::elements() const -> const std::vector<Element>&
{
    return pimpl->elements;
}

auto Phase::species() const -> const std::vector<Species>&
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
    return physicalState() == PhasePhysicalState::Solid;
}

auto Phase::standardThermoModel() const -> const PhaseStandardThermoModelFn&
{
    return pimpl->standard_thermo_model_fn;
}

auto Phase::activityModel() const -> const PhaseActivityModelFn&
{
    return pimpl->activity_model_fn;
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

auto Phase::standardThermoProps(double T, double P) const -> PhaseStandardThermoProps
{
    return pimpl->standard_thermo_model_fn(T, P);
}

auto Phase::activityProps(double T, double P, VectorConstRef n) const -> PhaseActivityProps
{
    return pimpl->activity_model_fn(T, P, n);
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
