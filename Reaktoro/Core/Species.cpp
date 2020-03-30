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

#include "Species.hpp"

// Reaktoro includes
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/Element.hpp>

namespace Reaktoro {

struct Species::Impl
{
    /// The name of the species such as `H2O(aq)`, `O2(g)`, `H+(aq)`.
    std::string name;

    /// The chemical formula of the species such as `H2O`, `O2`, `H+`.
    ChemicalFormula formula;

    /// The aggregate state of the species such as `aqueous`, `gaseous`, `liquid`, `solid`, etc..
    AggregateState aggregate_state;

    /// The tags of the species such as `organic`, `mineral`.
    std::vector<std::string> tags;

    /// The attached data whose type is known at runtime only and their ids.
    std::unordered_map<std::string, std::any> attached_data;

    /// Construct a default Species::Impl instance
    Impl()
    {}

    /// Construct a Species::Impl instance
    Impl(const ChemicalFormula& formula)
    : name(formula.str()), formula(formula)
    {}
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(const ChemicalFormula& formula)
: pimpl(new Impl(formula))
{}

auto Species::withName(std::string name) -> Species
{
    Species copy = clone();
    copy.pimpl->name = std::move(name);
    return copy;
}

auto Species::withFormula(const ChemicalFormula& formula) -> Species
{
    Species copy = clone();
    copy.pimpl->formula = formula;
    return copy;
}

auto Species::withAggregateState(AggregateState option) -> Species
{
    Species copy = clone();
    copy.pimpl->aggregate_state = option;
    return copy;
}

auto Species::withTags(std::vector<std::string> tags) -> Species
{
    Species copy = clone();
    copy.pimpl->tags = std::move(tags);
    return copy;
}

auto Species::withAttachedData(std::string id, std::any data) -> Species
{
    Species copy = clone();
    copy.pimpl->attached_data[id] = std::move(data);
    return copy;
}

auto Species::name() const -> std::string
{
    if(pimpl->name.empty())
        return formula().str();
    return pimpl->name;
}

auto Species::formula() const -> const ChemicalFormula&
{
    return pimpl->formula;
}

auto Species::charge() const -> double
{
    return formula().charge();
}

auto Species::molarMass() const -> double
{
    return formula().molarMass();
}

auto Species::aggregateState() const -> AggregateState
{
    if(pimpl->aggregate_state == AggregateState::Undefined)
        return identifyAggregateState(name());
    return pimpl->aggregate_state;
}

auto Species::elements() const -> const std::vector<std::pair<Element, double>>&
{
    return formula().elements();
}

auto Species::elementCoefficient(const std::string& symbol) const -> double
{
    return formula().coefficient(symbol);
}

auto Species::tags() const -> const std::vector<std::string>&
{
    return pimpl->tags;
}

auto Species::attachedData(std::string id) const -> std::optional<std::any>
{
    const auto iter = pimpl->attached_data.find(id);
    if(iter != pimpl->attached_data.end()) return iter->second;
    return {};
}

auto Species::attachedData() const -> const std::unordered_map<std::string, std::any>&
{
    return pimpl->attached_data;
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
