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

#include "Species.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>

namespace Reaktoro {
namespace detail {

/// Return the name of the species from its formula by removing suffix such as (aq), (s, calcite), etc..
auto removeSuffix(const String& formula) -> String
{
    return splitSpeciesNameSuffix(formula).first; // remove suffix
}

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

    /// The aggregate state of the species.
    AggregateState aggregatestate;

    /// The tags of the species such as `organic`, `mineral`.
    Strings tags;

    /// The attached data whose type is known at runtime only.
    std::any attacheddata;

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
      aggregatestate(identifyAggregateState(formula))
    {}
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(const String& formula)
: pimpl(new Impl(formula))
{}

auto Species::withName(const String& name) -> Species
{
    Species copy = clone();
    copy.pimpl->name = name;
    return copy;
}

auto Species::withFormula(const String& formula) -> Species
{
    Species copy = clone();
    copy.pimpl->formula = formula;
    return copy;
}

auto Species::withSubstance(const String& substance) -> Species
{
    Species copy = clone();
    copy.pimpl->substance = substance;
    return copy;
}

auto Species::withElements(const ElementalComposition& elements) -> Species
{
    Species copy = clone();
    copy.pimpl->elements = elements;
    return copy;
}

auto Species::withCharge(double charge) -> Species
{
    Species copy = clone();
    copy.pimpl->charge = charge;
    return copy;
}

auto Species::withAggregateState(AggregateState option) -> Species
{
    Species copy = clone();
    copy.pimpl->aggregatestate = option;
    return copy;
}

auto Species::withTags(const Strings& tags) -> Species
{
    Species copy = clone();
    copy.pimpl->tags = tags;
    return copy;
}

auto Species::withAttachedData(std::any data) -> Species
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

auto Species::charge() const -> double
{
    return pimpl->charge;
}

auto Species::molarMass() const -> double
{
    return elements().molarMass();
}

auto Species::aggregateState() const -> AggregateState
{
    if(pimpl->aggregatestate == AggregateState::Undefined)
        return identifyAggregateState(name());
    return pimpl->aggregatestate;
}

auto Species::elements() const -> const ElementalComposition&
{
    return pimpl->elements;
}

auto Species::tags() const -> const Strings&
{
    return pimpl->tags;
}

auto Species::attachedData() const -> const std::any&
{
    return pimpl->attacheddata;
}

auto Species::clone() const -> Species
{
    Species species;
    *species.pimpl = *pimpl;
    return species;
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
