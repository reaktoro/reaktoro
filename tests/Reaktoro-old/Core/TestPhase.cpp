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

#include "TestPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

namespace Reaktoro {
namespace {

auto test_Phase() -> void
{
	std::vector<Species> species(2);
	species[0].setName("A");
	species[1].setName("B");

	Phase phase;
	phase.setName("Aqueous");
	phase.setSpecies(species);

    ASSERT_EQUAL("Aqueous", phase.name());
    ASSERT_EQUAL(species.size(), phase.species().size());
    ASSERT_EQUAL(species[0].name(), phase.species()[0].name());
    ASSERT_EQUAL(species[1].name(), phase.species()[1].name());
}

auto test_numSpecies() -> void
{
	Phase phase;
	phase.setSpecies(std::vector<Species>(3));
	ASSERT_EQUAL(3, numSpecies(phase));
}

auto test_indexSpecies() -> void
{
	std::vector<Species> species(3);
	species[0].setName("A");
	species[1].setName("B");
	species[2].setName("C");

	Phase phase;
	phase.setSpecies(species);

	ASSERT_EQUAL(0, indexSpecies(phase, "A"));
	ASSERT_EQUAL(1, indexSpecies(phase, "B"));
	ASSERT_EQUAL(2, indexSpecies(phase, "C"));
	ASSERT_EQUAL(numSpecies(phase), indexSpecies(phase, "H"));
}

auto test_containsSpecies() -> void
{
	std::vector<Species> species(3);
	species[0].setName("A");
	species[1].setName("B");
	species[2].setName("C");

	Phase phase;
	phase.setSpecies(species);

	ASSERT(containsSpecies(phase, "A"));
	ASSERT(containsSpecies(phase, "B"));
	ASSERT(containsSpecies(phase, "C"));
	ASSERT(not containsSpecies(phase, "H"));
}

auto test_phaseNames() -> void
{
	std::vector<Phase> phases(2);
	phases[0].setName("Phase0");
	phases[1].setName("Phase1");
	const std::vector<std::string> phase_names = {"Phase0", "Phase1"};
	ASSERT_EQUAL(phase_names, phaseNames(phases));
}

auto test_phaseSpeciesInterpolatedThermoProperties() -> void
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
    thermo_model.heat_capacity = thermo_property_fn;

    std::vector<Species> species(2);
    species[0].setThermoModel(thermo_model);
    species[1].setThermoModel(thermo_model);

    Phase phase;
    phase.setSpecies(species);

    ASSERT_EQUAL(thermo_properties, enthalpies(phase, 300, 1));
    ASSERT_EQUAL(thermo_properties, entropies(phase, 300, 1));
    ASSERT_EQUAL(thermo_properties, standardGibbsEnergies(phase, 300, 1));
    ASSERT_EQUAL(thermo_properties, heatCapacitiesCp(phase, 300, 1));
    ASSERT_EQUAL(thermo_properties, helmholtzEnergies(phase, 300, 1));
    ASSERT_EQUAL(thermo_properties, internalEnergies(phase, 300, 1));
    ASSERT_EQUAL(thermo_properties, volumes(phase, 300, 1));
}

auto test_moleFractions() -> void
{
	Vector n = {2.0, 8.0};
	ChemicalVector x(2, 2);
	x.row(0) = ChemicalScalar(0.2, 0.0, 0.0, Vector{+0.08, -0.02});
	x.row(1) = ChemicalScalar(0.8, 0.0, 0.0, Vector{-0.08, +0.02});
    Phase phase;
    phase.setSpecies(std::vector<Species>(2));
    ChemicalVector x_actual = moleFractions(phase, n);
    const double eps = 1.e-16;
	ASSERT_EQUAL_DELTA(x.val[0], x_actual.val[0], eps);
	ASSERT_EQUAL_DELTA(x.val[1], x_actual.val[1], eps);
	ASSERT_EQUAL_DELTA(x.ddn[0], x_actual.ddn[0], eps);
	ASSERT_EQUAL_DELTA(x.ddn[1], x_actual.ddn[1], eps);
	ASSERT_EQUAL_DELTA(x.ddn[2], x_actual.ddn[2], eps);
	ASSERT_EQUAL_DELTA(x.ddn[3], x_actual.ddn[3], eps);
}

auto test_phaseThermoModels() -> void
{
	Vector n = {2.0, 8.0};
	ChemicalScalar rho(1000.0, -100.0, +200.0, Vector{0.2, 0.3});
	ChemicalVector c(2, 2);
	c.row(0) = ChemicalScalar(0.2, 0.0, 0.0, Vector{+0.08, -0.02});
	c.row(1) = ChemicalScalar(0.8, 0.0, 0.0, Vector{-0.08, +0.02});
	PhaseThermoModel thermo_model;
	thermo_model.concentration = [=](const auto&) { return c; };
	thermo_model.activity = [=](Temperature, Pressure, const auto&) { return c; };
	thermo_model.density = [=](Temperature, Pressure, const auto&) { return rho; };
	Phase phase;
    phase.setSpecies(std::vector<Species>(2));
    phase.setThermoModel(thermo_model);
	ASSERT_EQUAL(c, concentrations(phase, n));
	ASSERT_EQUAL(c, activities(phase, 300.0, 1.0, n));
	ASSERT_EQUAL(rho, density(phase, 300.0, 1.0, n));
}

} // namespace

auto testSuitePhase() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_Phase);
	s += CUTE(test_numSpecies);
	s += CUTE(test_indexSpecies);
	s += CUTE(test_containsSpecies);
	s += CUTE(test_phaseNames);
	s += CUTE(test_phaseSpeciesInterpolatedThermoProperties);
	s += CUTE(test_molarFractions);
	s += CUTE(test_phaseThermoModels);

    return s;
}

} // namespace Reaktoro
