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

// C++ includes
#include <cctype>
#include <set>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace {

/// Return the molar mass a species given its elemental composition.
auto molarMassFromElements(const std::map<Element, double>& elements) -> double
{
    double res = 0.0;
    for(auto pair : elements)
        res += pair.first.molarMass() * pair.second;
    return res;
}

/// Return the electrical charge of a species given its elemental composition.
auto chargeFromElements(const std::map<Element, double>& elements) -> double
{
    for(auto pair : elements)
        if(pair.first.name() == "Z")
            return pair.second;
    return 0.0;
}

} // namespace

struct Species::Impl
{
    /// The name of the chemical species
    std::string name;

    /// The chemical formula of the chemical species
    std::string formula;

    /// The elements that compose the chemical species and their coefficients
    std::map<Element, double> elements;

    /// The molar mass of the chemical species (in units of kg/mol)
    double molar_mass;

    /// The electrical charge of the chemical species
    double charge = 0.0;
};

Species::Species()
: pimpl(new Impl())
{}

auto Species::setName(std::string name) -> void
{
    pimpl->name = name;
}

auto Species::setFormula(std::string formula) -> void
{
    pimpl->formula = formula;
}

auto Species::setElements(const std::map<Element, double>& elements) -> void
{
    pimpl->elements = elements;
    pimpl->molar_mass = molarMassFromElements(elements);
    pimpl->charge = chargeFromElements(elements);
}

auto Species::numElements() const -> unsigned
{
    return elements().size();
}

auto Species::name() const -> std::string
{
    return pimpl->name;
}

auto Species::formula() const -> std::string
{
    return pimpl->formula;
}

auto Species::elements() const -> const std::map<Element, double>&
{
    return pimpl->elements;
}

auto Species::molarMass() const -> double
{
    return pimpl->molar_mass;
}

auto Species::charge() const -> double
{
    return pimpl->charge;
}

auto Species::elementCoefficient(std::string element) const -> double
{
    for(const auto& pair : elements())
        if(element == pair.first.name())
            return pair.second;
    return 0.0;
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
