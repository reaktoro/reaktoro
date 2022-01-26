// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "Material.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>

namespace Reaktoro {

struct Material::Impl
{
    /// The chemical formula of the substance (e.g., `HCO3-`).
    String formula;

    /// The element symbols and their coefficients (e.g., `{{"H", 1}, {"C", 1}, {"O", 3}}` for `HCO3-`).
    Pairs<String, double> elements;

    /// The electric charge in the chemical formula (e.g., `-1` for `HCO3-`).
    double charge = {};

    /// Construct an object of type Impl.
    Impl()
    {}

    /// Construct an object of type Impl with given formula.
    Impl(String formula)
    : formula(formula), elements(parseChemicalFormula(formula)), charge(parseElectricCharge(formula))
    {}

    /// Construct an object of type Impl with given data.
    Impl(String formula, Pairs<String, double> elements, double charge)
    : formula(formula), elements(elements), charge(charge)
    {}

    /// Return the symbols of the elements.
    auto symbols() const -> Strings
    {
        return vectorize(elements, RKT_LAMBDA(pair, pair.first));
    }

    /// Return the coefficients of the elements.
    auto coefficients() const -> Vec<double>
    {
        return vectorize(elements, RKT_LAMBDA(pair, pair.second));
    }

    /// Return the coefficient of an element symbol in the chemical formula.
    auto coefficient(const String& symbol) const -> double
    {
        for(auto const& element : elements)
            if(element.first == symbol)
                return element.second;
        return 0.0;
    }
};

Material::Material()
: pimpl(new Impl())
{}

Material::Material(const char* formula)
: Material(String(formula))
{}

Material::Material(String formula)
: pimpl(new Impl(formula))
{}

Material::Material(String formula, Pairs<String, double> symbols, double charge)
: pimpl(new Impl(formula, symbols, charge))
{}

auto Material::str() const -> const String&
{
    return pimpl->formula;
}

auto Material::elements() const -> const Pairs<String, double>&
{
    return pimpl->elements;
}

auto Material::symbols() const -> Strings
{
    return pimpl->symbols();
}

auto Material::coefficients() const -> Vec<double>
{
    return pimpl->coefficients();
}

auto Material::coefficient(const String& symbol) const -> double
{
    return pimpl->coefficient(symbol);
}

auto Material::charge() const -> double
{
    return pimpl->charge;
}

auto Material::equivalent(const Material& other) const -> bool
{
    return elements().size() == other.elements().size() &&
        contained(elements(), other.elements()) &&
        charge() == other.charge();
}

auto Material::equivalent(const Material& f1, const Material& f2) -> bool
{
    return f1.equivalent(f2);
}

Material::operator String() const
{
    return str();
}

Material::operator Pairs<String, double>() const
{
    return pimpl->elements;
}

auto operator<(const Material& lhs, const Material& rhs) -> bool
{
    return lhs.str() < rhs.str();
}

auto operator==(const Material& lhs, const Material& rhs) -> bool
{
    return lhs.str() == rhs.str();
}

} // namespace Reaktoro
