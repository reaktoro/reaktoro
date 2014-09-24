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

#include "SpeciesUtils.hpp"

// C++ includes
#include <algorithm>

// Reaktor includes
#include <Reaktor/Common/ThermoProperty.hpp>
#include <Reaktor/Common/ThermoProperties.hpp>
#include <Reaktor/Core/Species.hpp>

namespace Reaktor {

auto numElements(const Species& species) -> unsigned
{
	return species.elements().size();
}

auto containsElement(const Species& species, const std::string& element) -> bool
{
    return species.elements().count(element);
}

auto elementNames(const Species& species) -> std::vector<std::string>
{
    std::vector<std::string> names;
    names.reserve(species.elements().size());
    for(const auto& pair : species.elements())
        names.push_back(pair.first);
    return names;
}

auto elementAtoms(const Species& species) -> std::vector<double>
{
    std::vector<double> coefficients;
    coefficients.reserve(species.elements().size());
    for(const auto& pair : species.elements())
        coefficients.push_back(pair.second);
    return coefficients;
}

auto elementAtoms(const Species& species, const std::string& element) -> double
{
    const auto& iter = species.elements().find(element);
    return iter != species.elements().end() ? iter->second : 0.0;
}

template<typename PropertyFunction>
auto properties(const std::vector<Species>& species, double T, double P, PropertyFunction func) -> ThermoProperties
{
    const unsigned nspecies = species.size();
    Vector val(nspecies), ddt(nspecies), ddp(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
    {
        ThermoProperty prop = func(species[i], T, P);
        val[i] = prop.val();
        ddt[i] = prop.ddt();
        ddp[i] = prop.ddp();
    }
    return ThermoProperties(std::move(val), std::move(ddt), std::move(ddp));
}

auto volume(const Species& species, double T, double P) -> ThermoProperty
{
    return species.thermoModel().volume(T, P);
}

auto volumes(const std::vector<Species>& species, double T, double P) -> ThermoProperties
{
    return properties(species, T, P, volume);
}

auto entropy(const Species& species, double T, double P) -> ThermoProperty
{
    return species.thermoModel().entropy(T, P);
}

auto entropies(const std::vector<Species>& species, double T, double P) -> ThermoProperties
{
    return properties(species, T, P, entropy);
}

auto helmholtzEnergy(const Species& species, double T, double P) -> ThermoProperty
{
    return species.thermoModel().helmholtz_energy(T, P);
}

auto helmholtzEnergies(const std::vector<Species>& species, double T, double P) -> ThermoProperties
{
    return properties(species, T, P, helmholtzEnergy);
}

auto internalEnergy(const Species& species, double T, double P) -> ThermoProperty
{
    return species.thermoModel().internal_energy(T, P);
}

auto internalEnergies(const std::vector<Species>& species, double T, double P) -> ThermoProperties
{
    return properties(species, T, P, internalEnergy);
}

auto enthalpy(const Species& species, double T, double P) -> ThermoProperty
{
    return species.thermoModel().enthalpy(T, P);
}

auto enthalpies(const std::vector<Species>& species, double T, double P) -> ThermoProperties
{
    return properties(species, T, P, enthalpy);
}

auto gibbsEnergy(const Species& species, double T, double P) -> ThermoProperty
{
    return species.thermoModel().gibbs_energy(T, P);
}

auto gibbsEnergies(const std::vector<Species>& species, double T, double P) -> ThermoProperties
{
    return properties(species, T, P, gibbsEnergy);
}

auto heatCapacityCp(const Species& species, double T, double P) -> ThermoProperty
{
    return species.thermoModel().heat_capacity_cp(T, P);
}

auto heatCapacitiesCp(const std::vector<Species>& species, double T, double P) -> ThermoProperties
{
    return properties(species, T, P, heatCapacityCp);
}

auto names(const std::vector<Species>& species) -> std::vector<std::string>
{
    std::vector<std::string> names(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        names[i] = species[i].name();
    return names;
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
