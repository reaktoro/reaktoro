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

// Cute includes
#include <cute/cute.h>
#include <cute/ide_listener.h>
#include <cute/cute_runner.h>

// Reaktor includes
#include <Reaktor/Common/ThermoProperties.hpp>
#include <Reaktor/Common/ThermoProperty.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/SpeciesUtils.hpp>
using namespace Reaktor;

#define ASSERT_EQUAL_ARMA(expected,actual) ASSERT(arma::all(expected==actual))

auto elements      = std::map<std::string, double>{{"A", 1}, {"B", 2}, {"C", 3}};
auto element_names = std::vector<std::string>{"A", "B", "C"};
auto element_atoms = std::vector<double>{1, 2, 3};

auto createSpeciesThermoModel() -> SpeciesThermoModel
{
    SpeciesThermoModel thermo_model;
    thermo_model.enthalpy         = [](double, double) { return ThermoProperty(1.0, 0.0, 0.0); };
    thermo_model.entropy          = [](double, double) { return ThermoProperty(2.0, 0.0, 0.0); };
    thermo_model.gibbs_energy     = [](double, double) { return ThermoProperty(3.0, 0.0, 0.0); };
    thermo_model.heat_capacity_cp = [](double, double) { return ThermoProperty(4.0, 0.0, 0.0); };
    thermo_model.helmholtz_energy = [](double, double) { return ThermoProperty(5.0, 0.0, 0.0); };
    thermo_model.internal_energy  = [](double, double) { return ThermoProperty(6.0, 0.0, 0.0); };
    thermo_model.volume           = [](double, double) { return ThermoProperty(7.0, 0.0, 0.0); };
    return thermo_model;
}

auto createSpecies() -> Species
{
    Species species;
    species.setName("AB2C3-");
    species.setCharge(-1.0);
    species.setElements(elements);
    species.setFormula("AB2C3");
    species.setMolarMass(100.0);
    species.setThermoModel(createSpeciesThermoModel());
    return species;
}

auto createSpeciesSet() -> std::vector<Species>
{
    std::vector<Species> species_set(2);
    species_set[0].setName("A");
    species_set[1].setName("B");
    species_set[0].setCharge(2.0);
    species_set[1].setCharge(3.0);
    species_set[0].setMolarMass(20.0);
    species_set[1].setMolarMass(30.0);
    species_set[0].setThermoModel(createSpeciesThermoModel());
    species_set[1].setThermoModel(createSpeciesThermoModel());
    return species_set;
}

auto species               = createSpecies();
auto species_set           = createSpeciesSet();
auto species_names         = std::vector<std::string>{"A", "B"};
auto species_charges       = Vector{2.0, 3.0};
auto species_molar_massses = Vector{20.0, 30.0};
auto enthalpy_vals         = Vector{1.0, 1.0};
auto entropy_vals          = Vector{2.0, 2.0};
auto gibbs_energy_vals     = Vector{3.0, 3.0};
auto heat_capacity_cp_vals = Vector{4.0, 4.0};
auto helmholtz_energy_vals = Vector{5.0, 5.0};
auto internal_energy_vals  = Vector{6.0, 6.0};
auto volume_vals           = Vector{7.0, 7.0};

auto test_numElements() -> void
{
    ASSERT_EQUAL(3, numElements(species));
}

auto test_containsElement() -> void
{
    ASSERT(containsElement(species, "A"));
    ASSERT(containsElement(species, "B"));
    ASSERT(containsElement(species, "C"));
    ASSERT(not containsElement(species, "H"));
}

auto test_elementNames() -> void
{
    ASSERT_EQUAL(element_names, elementNames(species));
}

auto test_elementAtoms() -> void
{
    ASSERT_EQUAL(element_atoms, elementAtoms(species));
    ASSERT_EQUAL(1, elementAtoms(species, "A"));
    ASSERT_EQUAL(2, elementAtoms(species, "B"));
    ASSERT_EQUAL(3, elementAtoms(species, "C"));
    ASSERT_EQUAL(0, elementAtoms(species, "H"));
}

auto test_enthalpy() -> void
{
    ASSERT_EQUAL(1.0, enthalpy(species, 300, 1).val());
}

auto test_entropy() -> void
{
    ASSERT_EQUAL(2.0, entropy(species, 300, 1).val());
}

auto test_gibbsEnergy() -> void
{
    ASSERT_EQUAL(3.0, gibbsEnergy(species, 300, 1).val());
}

auto test_heatCapacityCp() -> void
{
    ASSERT_EQUAL(4.0, heatCapacityCp(species, 300, 1).val());
}

auto test_helmholtzEnergy() -> void
{
    ASSERT_EQUAL(5.0, helmholtzEnergy(species, 300, 1).val());
}

auto test_internalEnergy() -> void
{
    ASSERT_EQUAL(6.0, internalEnergy(species, 300, 1).val());
}

auto test_volume() -> void
{
    ASSERT_EQUAL(7.0, volume(species, 300, 1).val());
}

auto test_speciesNames() -> void
{
    ASSERT_EQUAL(species_names, speciesNames(species_set));    
}

auto test_speciesCharges() -> void
{
    ASSERT_EQUAL_ARMA(species_charges, speciesCharges(species_set));    
}

auto test_speciesMolarMasses() -> void
{
    ASSERT_EQUAL_ARMA(species_molar_massses, speciesMolarMasses(species_set));    
}

auto test_enthalpies() -> void
{
    ASSERT_EQUAL_ARMA(enthalpy_vals, enthalpies(species_set, 300, 1).val());    
}

auto test_entropies() -> void
{
    ASSERT_EQUAL_ARMA(entropy_vals, entropies(species_set, 300, 1).val());    
}

auto test_gibbsEnergies() -> void
{
    ASSERT_EQUAL_ARMA(gibbs_energy_vals, gibbsEnergies(species_set, 300, 1).val());    
}

auto test_heatCapacitiesCp() -> void
{
    ASSERT_EQUAL_ARMA(heat_capacity_cp_vals, heatCapacitiesCp(species_set, 300, 1).val());    
}

auto test_helmholtzEnergies() -> void
{
    ASSERT_EQUAL_ARMA(helmholtz_energy_vals, helmholtzEnergies(species_set, 300, 1).val());    
}

auto test_internalEnergies() -> void
{
    ASSERT_EQUAL_ARMA(internal_energy_vals, internalEnergies(species_set, 300, 1).val());    
}

auto test_volumes() -> void
{
    ASSERT_EQUAL_ARMA(volume_vals, volumes(species_set, 300, 1).val());    
}

int main(int argc, char **argv)
{
    cute::suite s;

    s.push_back(CUTE(test_numElements));
    s.push_back(CUTE(test_containsElement));
    s.push_back(CUTE(test_elementNames));
    s.push_back(CUTE(test_elementAtoms));
    s.push_back(CUTE(test_enthalpy));
    s.push_back(CUTE(test_entropy));
    s.push_back(CUTE(test_gibbsEnergy));
    s.push_back(CUTE(test_heatCapacityCp));
    s.push_back(CUTE(test_helmholtzEnergy));
    s.push_back(CUTE(test_internalEnergy));
    s.push_back(CUTE(test_volume));
    s.push_back(CUTE(test_speciesNames));
    s.push_back(CUTE(test_speciesCharges));
    s.push_back(CUTE(test_speciesMolarMasses));
    s.push_back(CUTE(test_enthalpies));
    s.push_back(CUTE(test_entropies));
    s.push_back(CUTE(test_gibbsEnergies));
    s.push_back(CUTE(test_heatCapacitiesCp));
    s.push_back(CUTE(test_helmholtzEnergies));
    s.push_back(CUTE(test_internalEnergies));
    s.push_back(CUTE(test_volumes));

    cute::ide_listener<> lis;
    cute::makeRunner(lis)(s, "SpeciesUtils");
}
