// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "GeneralSpecies.hpp"

// Reaktoro includes
#include <Reaktoro/Core/Element.hpp>

namespace Reaktoro {

struct GeneralSpecies::Impl
{
    /// The name of the chemical species
    std::string name;

    /// The chemical formula of the chemical species
    std::string formula;

    /// The elements that compose the chemical species and their coefficients
    std::map<Element, double> elements;

    /// The molar mass of the chemical species (in units of kg/mol)
    double molar_mass;
};

GeneralSpecies::GeneralSpecies()
: pimpl(new Impl())
{}

GeneralSpecies::GeneralSpecies(const GeneralSpecies& other)
: pimpl(new Impl(*other.pimpl))
{}

GeneralSpecies::~GeneralSpecies()
{}

auto GeneralSpecies::operator=(GeneralSpecies other) -> GeneralSpecies&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto GeneralSpecies::setName(std::string name) -> void
{
    pimpl->name = name;
}

auto GeneralSpecies::setFormula(std::string formula) -> void
{
    pimpl->formula = formula;
}

auto GeneralSpecies::setElements(const std::map<Element, double>& elements) -> void
{
    pimpl->elements = elements;
}

auto GeneralSpecies::setMolarMass(double value) -> void
{
    pimpl->molar_mass = value;
}

auto GeneralSpecies::numElements() const -> unsigned
{
    return elements().size();
}

auto GeneralSpecies::name() const -> const std::string&
{
    return pimpl->name;
}

auto GeneralSpecies::formula() const -> const std::string&
{
    return pimpl->formula;
}

auto GeneralSpecies::elements() const -> const std::map<Element, double>&
{
    return pimpl->elements;
}

auto GeneralSpecies::molarMass() const -> double
{
    return pimpl->molar_mass;
}

auto GeneralSpecies::elementCoefficient(std::string element) const -> double
{
    for(const auto& pair : elements())
        if(element == pair.first.name())
            return pair.second;
    return 0.0;
}

auto operator<(const GeneralSpecies& lhs, const GeneralSpecies& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const GeneralSpecies& lhs, const GeneralSpecies& rhs) -> bool
{

    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
