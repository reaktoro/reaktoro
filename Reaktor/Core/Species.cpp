// Reaktor is a C++ library for computational reaction modelling.
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

#include <Reaktor/Core/CoreUtils.hpp>
#include "Species.hpp"

// C++ includes
#include <set>

// Reaktor includes

namespace Reaktor {

struct Species::Impl
{
    SpeciesData data;
};

Species::Species()
: pimpl(new Impl())
{}

Species::Species(const SpeciesData& data)
: pimpl(new Impl{data})
{}

Species::Species(std::string name, std::string formula, std::vector<Element> elements, std::vector<double> atoms, double charge, double molar_mass)
: Species(SpeciesData{name, formula, elements, atoms, charge, molar_mass})
{}

auto Species::numElements() const -> unsigned
{
    return elements().size();
}

auto Species::name() const -> const std::string&
{
    return pimpl->data.name;
}

auto Species::formula() const -> const std::string&
{
    return pimpl->data.formula;
}

auto Species::elements() const -> const std::vector<Element>&
{
    return pimpl->data.elements;
}

auto Species::atoms() const -> const std::vector<double>&
{
    return pimpl->data.atoms;
}

auto Species::charge() const -> double
{
    return pimpl->data.charge;
}

auto Species::molarMass() const -> double
{
    return pimpl->data.molar_mass;
}

auto Species::standardGibbsEnergy(double T, double P) const -> ThermoScalar
{
    return pimpl->data.standard_gibbs_energy(T, P);
}

auto Species::standardHelmholtzEnergy(double T, double P) const -> ThermoScalar
{
    return pimpl->data.standard_helmholtz_energy(T, P);
}

auto Species::standardInternalEnergy(double T, double P) const -> ThermoScalar
{
    return pimpl->data.standard_internal_energy(T, P);
}

auto Species::standardEnthalpy(double T, double P) const -> ThermoScalar
{
    return pimpl->data.standard_enthalpy(T, P);
}

auto Species::standardEntropy(double T, double P) const -> ThermoScalar
{
    return pimpl->data.standard_entropy(T, P);
}

auto Species::standardVolume(double T, double P) const -> ThermoScalar
{
    return pimpl->data.standard_volume(T, P);
}

auto Species::standardHeatCapacity(double T, double P) const -> ThermoScalar
{
    return pimpl->data.standard_heat_capacity(T, P);
}

auto operator<(const Species& lhs, const Species& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Species& lhs, const Species& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

auto atoms(const Element& element, const Species& species) -> double
{
    Index idx = index(element, species.elements());
    return idx < species.elements().size() ? species.atoms()[idx] : 0.0;
}

auto collectElements(const std::vector<Species>& species) -> std::vector<Element>
{
    std::set<Element> elements;
    for(const Species& iter : species)
        elements.insert(iter.elements().begin(), iter.elements().end());
    return std::vector<Element>(elements.begin(), elements.end());
}

auto charges(const std::vector<Species>& species) -> Vector
{
    Vector charges(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        charges[i] = species[i].charge();
    return charges;
}

auto molarMasses(const std::vector<Species>& species) -> Vector
{
    Vector molar_masses(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        molar_masses[i] = species[i].molarMass();
    return molar_masses;
}

} // namespace Reaktor
