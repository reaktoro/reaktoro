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
    AggregateState aggregatestate;

    /// The base function that computes the standard thermodynamic properties of the species (if any).
    /// The result of this function may be overwritten if a formation reaction with thermo props are given.
    StandardThermoPropsFn standard_thermo_props_fn;

    /// The function used by method props(T, P) to compute the standard thermo props of the species.
    /// This function combines the function `standard_thermo_props_fn` with the
    /// functions `reaction.standardGibbsEnergy()` and `reaction.standardEnthalpy()`
    /// so that if thermodynamic data is provided in terms of formation reaction, namely
    /// lg(K) and delta(H0), these can be used together with others provided directly as standard
    /// thermodynamic properties (e.g., V0, Cp0, CV0).
    StandardThermoPropsFn props_fn;

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
      aggregatestate(identifyAggregateState(formula))
    {}

    /// Initialize the final function for standard thermodynamic property evaluations of the species.
    auto initializePropsFn()
    {
        const auto zeroprops = [](real T, real P) { return StandardThermoProps{}; }; // return zeros for all standard properties
        const auto baseprops = standard_thermo_props_fn ? standard_thermo_props_fn : zeroprops; // ensure non-empty standard thermo props fn
        if(reaction.reactants().size())
        {
            const auto G0 = reaction.standardGibbsEnergyFn();
            const auto H0 = reaction.standardEnthalpyFn();
            props_fn = [=](real T, real P)
            {
                StandardThermoProps props = baseprops(T, P);
                if(G0) props.G0 = G0(T, P); // overwrite G0 value if reaction has lgK data
                if(H0) props.H0 = H0(T, P); // overwrite H0 value if reaction has dH0 data
                return props;
            };
        }
        else props_fn = baseprops;

        // Use memoization to ensure fast re-evaluation when new
        // T and P inputs are identical to last ones.
        if(props_fn) props_fn = memoizeLast(props_fn);
    }

    /// Return the standard thermodynamic properties of the species at given temperature (in K) and pressure (in Pa).
    auto props(real T, real P) const -> StandardThermoProps
    {
        error(!props_fn, "Cannot compute the standard thermodynamic properties of Species object with name ", name, ". \n"
            "To fix this error, use one of the methods below in this Species object: \n"
            "    1) Species::withStandardThermoPropsFn\n"
            "    2) Species::withStandardGibbsEnergyFn\n"
            "    3) Species::withStandardGibbsEnergy\n"
            "    4) Species::withFormationReaction");
        return props_fn(T, P);
    }
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(String formula)
: pimpl(new Impl(formula))
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
    copy.pimpl->aggregatestate = option;
    return copy;
}

auto Species::withFormationReaction(FormationReaction reaction) const -> Species
{
    Species copy = clone();
    copy.pimpl->reaction = std::move(reaction);
    copy.pimpl->initializePropsFn();
    return copy;
}

auto Species::withStandardGibbsEnergy(real value) const -> Species
{
    return withStandardGibbsEnergyFn([=](real T, real P) { return value; });
}

auto Species::withStandardGibbsEnergyFn(const Fn<real(real,real)>& fn) const -> Species
{
    return withStandardThermoPropsFn([=](real T, real P)
    {
        StandardThermoProps props = {};
        props.G0 = fn(T, P);
        return props;
    });
}

auto Species::withStandardThermoPropsFn(const StandardThermoPropsFn& fn) const -> Species
{
    Species copy = clone();
    copy.pimpl->standard_thermo_props_fn = memoizeLast(fn);
    copy.pimpl->initializePropsFn();
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
    if(pimpl->aggregatestate == AggregateState::Undefined)
        return identifyAggregateState(name());
    return pimpl->aggregatestate;
}

auto Species::reaction() const -> const FormationReaction&
{
    return pimpl->reaction;
}

auto Species::standardThermoPropsFn() const -> const StandardThermoPropsFn&
{
    return pimpl->standard_thermo_props_fn;
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
    return pimpl->props(T, P);
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
