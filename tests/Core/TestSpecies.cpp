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

#include "TestSpecies.hpp"

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

namespace Reaktor {
namespace {

auto test_Species() -> void
{
    Species species;
    species.setName("AB2C3-");
    species.setCharge(-1.0);
    species.setElements({"A", "B", "C"});
    species.setElementAtoms({1, 2, 3});
    species.setFormula("AB2C3-");
    species.setMolarMass(100.0);

    ASSERT_EQUAL("AB2C3-" , species.name());
    ASSERT_EQUAL(-1.0     , species.charge());
    ASSERT_EQUAL(3        , species.elements().size());
    ASSERT_EQUAL(3        , species.elementAtoms().size());
    ASSERT_EQUAL("A"      , species.elements()[0]);
    ASSERT_EQUAL("B"      , species.elements()[1]);
    ASSERT_EQUAL("C"      , species.elements()[2]);
    ASSERT_EQUAL(1        , species.elementAtoms()[0]);
    ASSERT_EQUAL(2        , species.elementAtoms()[1]);
    ASSERT_EQUAL(3        , species.elementAtoms()[2]);
    ASSERT_EQUAL("AB2C3-" , species.formula());
    ASSERT_EQUAL(100.0    , species.molarMass());
}

auto test_numElements() -> void
{
    Species species;
    species.setElements({"H", "O"});
    species.setElementAtoms({1, 2});
    ASSERT_EQUAL(2, numElements(species));
}

auto test_containsElement() -> void
{
    Species species;
    species.setElements({"H", "O"});
    ASSERT(containsElement(species, "H"));
    ASSERT(containsElement(species, "O"));
    ASSERT(not containsElement(species, "N"));
}

auto test_elementIndex() -> void
{
    Species species;
    species.setElements({"H", "O"});
    ASSERT_EQUAL(0, elementIndex(species, "H"));
    ASSERT_EQUAL(1, elementIndex(species, "O"));
    ASSERT_EQUAL(numElements(species), elementIndex(species, "N"));
}

auto test_speciesThermoVector() -> void
{
    ThermoScalar thermo_property(1.0, 2.0, 3.0);
    ThermoVector thermo_properties(Vector{1.0, 1.0}, Vector{2.0, 2.0}, Vector{3.0, 3.0});
    ThermoScalarFunction thermo_property_fn = [=](double,double) { return thermo_property; };
    SpeciesThermoModel thermo_model;
    thermo_model.gibbs_energy     = thermo_property_fn;
    thermo_model.helmholtz_energy = thermo_property_fn;
    thermo_model.internal_energy  = thermo_property_fn;
    thermo_model.enthalpy         = thermo_property_fn;
    thermo_model.entropy          = thermo_property_fn;
    thermo_model.volume           = thermo_property_fn;
    thermo_model.heat_capacity_cp = thermo_property_fn;

    std::vector<Species> species(2);
    species[0].setThermoModel(thermo_model);
    species[1].setThermoModel(thermo_model);

    ASSERT_EQUAL(thermo_property, enthalpy(species[0], 300, 1));
    ASSERT_EQUAL(thermo_property, entropy(species[0], 300, 1));
    ASSERT_EQUAL(thermo_property, gibbsEnergy(species[0], 300, 1));
    ASSERT_EQUAL(thermo_property, heatCapacityCp(species[0], 300, 1));
    ASSERT_EQUAL(thermo_property, helmholtzEnergy(species[0], 300, 1));
    ASSERT_EQUAL(thermo_property, internalEnergy(species[0], 300, 1));
    ASSERT_EQUAL(thermo_property, volume(species[0], 300, 1));

    ASSERT_EQUAL(thermo_properties, enthalpies(species, 300, 1));
    ASSERT_EQUAL(thermo_properties, entropies(species, 300, 1));
    ASSERT_EQUAL(thermo_properties, gibbsEnergies(species, 300, 1));
    ASSERT_EQUAL(thermo_properties, heatCapacitiesCp(species, 300, 1));
    ASSERT_EQUAL(thermo_properties, helmholtzEnergies(species, 300, 1));
    ASSERT_EQUAL(thermo_properties, internalEnergies(species, 300, 1));
    ASSERT_EQUAL(thermo_properties, volumes(species, 300, 1));
}

auto test_speciesNames() -> void
{
    std::vector<std::string> names = {"Species0", "Species1"};
    std::vector<Species> species(2);
    species[0].setName(names[0]);
    species[1].setName(names[1]);
    ASSERT_EQUAL(names, speciesNames(species));
}

auto test_speciesCharges() -> void
{
    Vector charges = {-1.0, +1.0};
    std::vector<Species> species(2);
    species[0].setCharge(charges[0]);
    species[1].setCharge(charges[1]);
    ASSERT(arma::all(charges == speciesCharges(species)));
}

auto test_speciesMolarMasses() -> void
{
    Vector molar_masses = {10.0, 20.0};
    std::vector<Species> species(2);
    species[0].setMolarMass(molar_masses[0]);
    species[1].setMolarMass(molar_masses[1]);
    ASSERT(arma::all(molar_masses == speciesMolarMasses(species)));
}

} // namespace

auto testSuiteSpecies() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_Species);
    s += CUTE(test_numElements);
    s += CUTE(test_containsElement);
    s += CUTE(test_elementIndex);
    s += CUTE(test_speciesThermoVector);
    s += CUTE(test_speciesNames);
    s += CUTE(test_speciesCharges);
    s += CUTE(test_speciesMolarMasses);

    return s;
}

} // namespace Reaktor
