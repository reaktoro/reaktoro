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
#include <Reaktoro/Singletons/PeriodicTable.hpp>

namespace Reaktoro {
namespace detail {

/// Create the Element objects in a Species object using PeriodicTable.
auto createElements(const ChemicalFormula& formula) -> Species::Elements
{
    Species::Elements elements;
    for(auto&& [symbol, coeff] : formula.symbols())
    {
        const auto element = PeriodicTable::elementWithSymbol(symbol);
        error(!element.has_value(), "Cannot construct Species object with formula",
            formula.str(), ". PeriodicTable contains no element with symbol ", symbol, ". "
            "Use method PeriodicTable::append (in C++) or PeriodicTable.append (in Python) "
            "to add a new Element with this symbol.");
        elements.emplace(element.value(), coeff);
    }
    return elements;
}

/// Create the Element objects in a Species object using PeriodicTable.
auto createElements(const Species::ElementSymbols& symbols) -> Species::Elements
{
    Species::Elements elements;
    for(auto&& [symbol, coeff] : symbols)
    {
        const auto element = PeriodicTable::elementWithSymbol(symbol);
        error(!element.has_value(), "Cannot proceed with Species::withElementSymbols. "
            "PeriodicTable contains no element with symbol ", symbol, ". "
            "Use method PeriodicTable::append (in C++) or PeriodicTable.append (in Python) "
            "to add a new Element with this symbol.");
        elements.emplace(element.value(), coeff);
    }
    return elements;
}

/// Return the element symbols and their coefficients in a Species::Elements container.
auto elementSymbols(const Species::Elements& elements)
{
    Species::ElementSymbols symbols;
    for(auto&& [element, coeff] : elements)
        symbols[element.symbol()] = coeff;
    return symbols;
}

/// Return the molar mass of a species.
auto molarMass(const Species::Elements& elements)
{
    double molar_mass = {};
    for(auto&& [element, coeff] : elements)
        molar_mass += element.molarMass() * coeff;
    return molar_mass;
}

/// Return the name of the species from its formula by removing suffix such as (aq), (s, calcite), etc..
auto removeSuffix(std::string formula) -> std::string
{
    return splitSpeciesNameSuffix(formula).first; // remove suffix
}

} // namespace detail

struct Species::Impl
{
    /// The name that uniquely identifies this species such as `H2O(aq)`, `O2(g)`, `H+(aq)`.
    std::string name;

    /// The chemical formula of the species such as `H2O`, `O2`, `H+`.
    std::string formula;

    /// The name of the underlying substance such as `H2O`, `WATER`, `CARBON-MONOXIDE`, `CO2`.
    std::string substance;

    /// The elements in the species and their coefficients.
    Elements elements;

    /// The electric charge of the species.
    double charge = {};

    /// The aggregate state of the species such as `aqueous`, `gaseous`, `liquid`, `solid`, etc..
    AggregateState aggregate_state;

    /// The tags of the species such as `organic`, `mineral`.
    std::vector<std::string> tags;

    /// The critical properties of the underlying substance of the species if available.
    std::optional<SubstanceCriticalProps> crprops;

    /// The attached data whose type is known at runtime only.
    std::any attached_data;

    /// Construct a default Species::Impl instance
    Impl()
    {}

    /// Construct a Species::Impl instance
    Impl(ChemicalFormula formula)
    : name(formula), formula(detail::removeSuffix(formula)), substance(detail::removeSuffix(formula)),
      elements(detail::createElements(formula)), charge(formula.charge())
    {
    }

    /// Return a ChemicalFormula object consistent with the elements and charge
    /// of the species, and not just based on the contents of string `formula`.
    auto createChemicalFormula() const -> ChemicalFormula
    {
        const auto symbols = detail::elementSymbols(elements);
        return ChemicalFormula(formula, symbols, charge);
    }
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(std::string formula)
: pimpl(new Impl(formula))
{}

auto Species::withName(std::string name) -> Species
{
    Species copy = clone();
    copy.pimpl->name = std::move(name);
    return copy;
}

auto Species::withFormula(std::string formula) -> Species
{
    Species copy = clone();
    copy.pimpl->formula = std::move(formula);
    return copy;
}

auto Species::withSubstance(std::string substance) -> Species
{
    Species copy = clone();
    copy.pimpl->substance = std::move(substance);
    return copy;
}

auto Species::withElements(Elements elements) -> Species
{
    Species copy = clone();
    copy.pimpl->elements = std::move(elements);
    return copy;
}

auto Species::withElementSymbols(ElementSymbols symbols) -> Species
{
    return withElements(detail::createElements(symbols));
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
    copy.pimpl->aggregate_state = option;
    return copy;
}

auto Species::withTags(std::vector<std::string> tags) -> Species
{
    Species copy = clone();
    copy.pimpl->tags = std::move(tags);
    return copy;
}

auto Species::withCriticalProps(const SubstanceCriticalProps& props) -> Species
{
    Species copy = clone();
    copy.pimpl->crprops = std::move(props);
    return copy;
}

auto Species::withAttachedData(std::any data) -> Species
{
    Species copy = clone();
    copy.pimpl->attached_data = std::move(data);
    return copy;
}

auto Species::name() const -> std::string
{
    if(pimpl->name.empty())
        return formula();
    return pimpl->name;
}

auto Species::formula() const -> ChemicalFormula
{
    return pimpl->createChemicalFormula();
}

auto Species::substance() const -> std::string
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
    return detail::molarMass(elements());
}

auto Species::aggregateState() const -> AggregateState
{
    if(pimpl->aggregate_state == AggregateState::Undefined)
        return identifyAggregateState(name());
    return pimpl->aggregate_state;
}

auto Species::elements() const -> const Elements&
{
    return pimpl->elements;
}

auto Species::tags() const -> const std::vector<std::string>&
{
    return pimpl->tags;
}

auto Species::criticalProps() const -> std::optional<SubstanceCriticalProps>
{
    if(pimpl->crprops)
        return pimpl->crprops;
    return CriticalProps::get({ substance(), formula(), name() });
}

auto Species::attachedData() const -> const std::any&
{
    return pimpl->attached_data;
}

auto Species::elementCoefficient(const std::string& symbol) const -> double
{
    for(auto&& [element, coeff] : elements())
        if(element.symbol() == symbol)
            return coeff;
    return 0.0;
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
