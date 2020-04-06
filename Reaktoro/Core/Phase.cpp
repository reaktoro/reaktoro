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
    StateOfMatter state = StateOfMatter::Solid;

    /// The list of Species instances defining the phase
    std::vector<Species> species;

    /// The list of Element instances in the phase
    std::vector<Element> elements;

    /// The standard thermodynamic model of the phase.
    StandardThermoModelFn standard_thermo_model_fn;

    /// The activity model of the phase.
    ActivityModelFn activity_model_fn;
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

auto Phase::setPhysicalState(StateOfMatter state) -> void
{
    pimpl->state = state;
}

auto Phase::setSpecies(const std::vector<Species>& species) -> void
{
    pimpl->species = species;
}

auto Phase::setStandardThermoModel(const StandardThermoModelFn& model) -> void
{
    pimpl->standard_thermo_model_fn = model;
}

auto Phase::setActivityModel(const ActivityModelFn& model) -> void
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

auto Phase::physicalState() const -> StateOfMatter
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
    return physicalState() == StateOfMatter::Solid;
}

auto Phase::standardThermoModel() const -> const StandardThermoModelFn&
{
    return pimpl->standard_thermo_model_fn;
}

auto Phase::activityModel() const -> const ActivityModelFn&
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

auto Phase::standardThermoProps(Index ispecies, real T, real P) const -> StandardThermoProps
{
    return pimpl->standard_thermo_model_fn(T, P, pimpl->species[ispecies]);
}

auto Phase::activityProps(real T, real P, VectorXrConstRef n) const -> ActivityProps
{
    ActivityProps props;
    props.ln_activity_coefficients.resize(numSpecies());
    props.ln_activities.resize(numSpecies());
    pimpl->activity_model_fn(props, T, P, n);
    return props;
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
