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
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

struct Phase::Impl
{
    /// The name of the phase
    std::string name;

    /// The physical state of the phase.
    StateOfMatter state = StateOfMatter::Solid;

    /// The list of Species instances defining the phase
    SpeciesList species;

    /// The standard thermodynamic model of the species in the phase.
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

auto Phase::withName(std::string name) -> Phase
{
    Phase copy = clone();
    copy.pimpl->name = std::move(name);
    return copy;
}

auto Phase::withSpecies(SpeciesList species) -> Phase
{
    Phase copy = clone();
    copy.pimpl->species = std::move(species);
    return copy;
}

auto Phase::withPhysicalState(StateOfMatter state) -> Phase
{
    Phase copy = clone();
    copy.pimpl->state = std::move(state);
    return copy;
}

auto Phase::withStandardThermoModel(StandardThermoModelFn model) -> Phase
{
    Phase copy = clone();
    copy.pimpl->standard_thermo_model_fn = std::move(model);
    return copy;
}

auto Phase::withActivityModel(ActivityModelFn model) -> Phase
{
    Phase copy = clone();
    copy.pimpl->activity_model_fn = std::move(model);
    return copy;
}

auto Phase::name() const -> std::string
{
    return pimpl->name;
}

auto Phase::physicalState() const -> StateOfMatter
{
    return pimpl->state;
}

auto Phase::species() const -> const SpeciesList&
{
    return pimpl->species;
}

auto Phase::species(Index idx) const -> const Species&
{
    return pimpl->species[idx];
}

auto Phase::standardThermoModel() const -> const StandardThermoModelFn&
{
    return pimpl->standard_thermo_model_fn;
}

auto Phase::activityModel() const -> const ActivityModelFn&
{
    return pimpl->activity_model_fn;
}

auto Phase::activityProps(real T, real P, VectorXrConstRef n) const -> ActivityProps
{
    assert(activityModel());
    ActivityProps props;
    props.ln_g.resize(species().size());
    props.ln_a.resize(species().size());
    activityModel()(props, T, P, n);
    return props;
}

auto Phase::standardThermoProps(real T, real P, Index ispecies) const -> StandardThermoProps
{
    assert(standardThermoModel());
    return standardThermoModel()(T, P, species(ispecies));
}

auto Phase::props(real T, real P, VectorXrConstRef n) const -> ChemicalPropsPhase
{
    ChemicalPropsPhase res;
    eval(res, T, P, n);
    return res;
}

auto Phase::eval(ChemicalPropsPhaseRef props, real T, real P, VectorXrConstRef n) const -> void
{

}

auto Phase::clone() const -> Phase
{
    Phase phase;
    *phase.pimpl = *pimpl;
    return phase;
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
