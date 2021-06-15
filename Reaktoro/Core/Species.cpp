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

#include "Species.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Memoization.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ElementalComposition.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>

namespace Reaktoro {
namespace detail {

/// Return the name of the species from its formula by removing suffix such as (aq), (s, calcite), etc..
auto removeSuffix(const String& formula) -> String
{
    return splitSpeciesNameSuffix(formula).first; // remove suffix
}

/// Return a default standard thermodynamic model for a species that always indicate failure to initialize the model.
auto defaultStandardThermoModel() -> StandardThermoModel
{
    return [](real T, real P) -> StandardThermoProps
    {
        errorif(true, "Cannot compute the standard thermodynamic properties of this Species object. \n"
            "To fix this error, use one of the methods below in this Species object: \n"
            "    1) Species::withStandardThermoModel\n"
            "    2) Species::withStandardGibbsEnergy\n"
            "    3) Species::withFormationReaction");
        return {};
    };
};

} // namespace detail

struct Species::Impl
{
    /// The name that uniquely identifies this species such as `H2O(aq)`, `O2(g)`, `H+(aq)`.
    String name;

    /// The chemical formula of the species such as `H2O`, `O2`, `H+`, `CO3--`, `CaMg(CO3)2`.
    String formula;

    /// The name of the underlying substance such as `H2O`, `WATER`, `CARBON-MONOXIDE`, `CO2`.
    String substance;

    /// The elemental composition of the species with its elements and respective coefficients.
    ElementalComposition elements;

    /// The electric charge of the species.
    double charge = 0.0;

    /// The formation reaction of the species and its thermodynamic properties (if any).
    FormationReaction reaction;

    /// The aggregate state of the species.
    AggregateState aggregate_state;

    /// The standard thermodynamic model function of the species (if any).
    StandardThermoModel propsfn = detail::defaultStandardThermoModel();

    /// The tags of the species such as `organic`, `mineral`.
    Strings tags;

    /// The attached data whose type is known at runtime only.
    Any attacheddata;

    /// Construct a default Species::Impl instance
    Impl()
    {}

    /// Construct a Species::Impl instance with given formula
    Impl(const ChemicalFormula& formula)
    : name(formula),
      formula(detail::removeSuffix(formula)),
      substance(detail::removeSuffix(formula)),
      elements(formula.elements()),
      charge(formula.charge()),
      aggregate_state(identifyAggregateState(formula))
    {}

    /// Construct a Species::Impl instance with given attributes
    Impl(const Attribs& attribs)
    {
        errorif(attribs.name.empty(), "Species::Attribs::name cannot be empty.");
        errorif(attribs.formula.empty(), "Species::Attribs::formula cannot be empty.");
        errorif(attribs.std_thermo_model && attribs.formation_reaction,
            "Species::Attribs for ", attribs.formula, " cannot contain both a "
            "FormationReaction object and a StandardThermoModel object.");
        name = attribs.name;
        formula = detail::removeSuffix(attribs.formula);
        substance = attribs.substance.value_or(formula);
        elements = attribs.elements.value_or(parseChemicalFormula(formula));
        charge = attribs.charge.value_or(parseElectricCharge(formula));
        aggregate_state = attribs.aggregate_state.value_or(identifyAggregateState(attribs.formula));
        tags = attribs.tags.value_or(Strings{});
        if(attribs.std_thermo_model)
        {
            propsfn = attribs.std_thermo_model.value();
            propsfn = propsfn.withMemoization();
        }
        if(attribs.formation_reaction)
        {
            reaction = attribs.formation_reaction.value();
            propsfn = reaction.createStandardThermoModel();
            propsfn = propsfn.withMemoization();
        }
    }
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(String formula)
: pimpl(new Impl(formula))
{}

Species::Species(const Attribs& attribs)
: pimpl(new Impl(attribs))
{}

auto Species::clone() const -> Species
{
    Species species;
    *species.pimpl = *pimpl;
    return species;
}

auto Species::withName(String name) const -> Species
{
    Species copy = clone();
    copy.pimpl->name = std::move(name);
    return copy;
}

auto Species::withFormula(String formula) const -> Species
{
    Species copy = clone();
    copy.pimpl->formula = std::move(formula);
    return copy;
}

auto Species::withSubstance(String substance) const -> Species
{
    Species copy = clone();
    copy.pimpl->substance = std::move(substance);
    return copy;
}

auto Species::withElements(ElementalComposition elements) const -> Species
{
    Species copy = clone();
    copy.pimpl->elements = std::move(elements);
    return copy;
}

auto Species::withCharge(double charge) const -> Species
{
    Species copy = clone();
    copy.pimpl->charge = charge;
    return copy;
}

auto Species::withAggregateState(AggregateState option) const -> Species
{
    Species copy = clone();
    copy.pimpl->aggregate_state = option;
    return copy;
}

auto Species::withFormationReaction(const FormationReaction& reaction) const -> Species
{
    Species copy = clone();
    copy.pimpl->reaction = reaction;
    copy = copy.withStandardThermoModel(reaction.createStandardThermoModel());
    return copy;
}

auto Species::withStandardGibbsEnergy(Param G0) const -> Species
{
    auto calcfn = [=](real T, real P)
    {
        StandardThermoProps props = {};
        props.G0 = G0;
        return props;
    };

    Params params = { G0 };

    return withStandardThermoModel({ calcfn, params });
}

auto Species::withStandardThermoModel(const StandardThermoModel& model) const -> Species
{
    Species copy = clone();
    copy.pimpl->propsfn = model.withMemoization();
    return copy;
}

auto Species::withTags(const Strings& tags) const -> Species
{
    Species copy = clone();
    copy.pimpl->tags = tags;
    return copy;
}

auto Species::withAttachedData(Any data) const -> Species
{
    Species copy = clone();
    copy.pimpl->attacheddata = data;
    return copy;
}

auto Species::name() const -> String
{
    if(pimpl->name.empty())
        return formula();
    return pimpl->name;
}

auto Species::formula() const -> ChemicalFormula
{
    return ChemicalFormula(pimpl->formula, pimpl->elements, pimpl->charge);
}

auto Species::substance() const -> String
{
    if(pimpl->substance.empty())
        return detail::removeSuffix(name());
    return pimpl->substance;
}

auto Species::elements() const -> const ElementalComposition&
{
    return pimpl->elements;
}

auto Species::charge() const -> double
{
    return pimpl->charge;
}

auto Species::aggregateState() const -> AggregateState
{
    if(pimpl->aggregate_state == AggregateState::Undefined)
        return identifyAggregateState(name());
    return pimpl->aggregate_state;
}

auto Species::reaction() const -> const FormationReaction&
{
    return pimpl->reaction;
}

auto Species::standardThermoModel() const -> const StandardThermoModel&
{
    return pimpl->propsfn;
}

auto Species::tags() const -> const Strings&
{
    return pimpl->tags;
}

auto Species::attachedData() const -> const Any&
{
    return pimpl->attacheddata;
}

auto Species::molarMass() const -> double
{
    return elements().molarMass();
}

auto Species::props(real T, real P) const -> StandardThermoProps
{
    return pimpl->propsfn(T, P);
}

auto operator<(const Species& lhs, const Species& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Species& lhs, const Species& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
