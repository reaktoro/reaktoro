// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// ReaktoroInterpreter includes
#include <ReaktoroInterpreterCpp/ParserUtils.hpp>
using namespace Reaktoro;

// cute includes
#include <cute/cute.h>
#include <cute/cute_runner.h>
#include <cute/ide_listener.h>

auto test_ValueUnitsParser() -> void
{
    const std::string s = "Temperature: 25 celsius";
    YAML::Node node = YAML::Load(s);
    ValueUnits x; node >> x;

    ASSERT_EQUAL(25, x.value);
    ASSERT_EQUAL("celsius", x.units);
}

auto test_EntityValueUnitsParser() -> void
{
    const std::string s = "SpeciesAmount: Calcite 100 g";
    YAML::Node node = YAML::Load(s);
    EntityValueUnits x; node >> x;

    ASSERT_EQUAL("Calcite", x.entity);
    ASSERT_EQUAL(100, x.value);
    ASSERT_EQUAL("g", x.units);
}

auto test_MixtureCompoundParser() -> void
{
    const std::string s = "1 kg H2O";
    YAML::Node node = YAML::Load(s);
    MixtureCompound x; node >> x;

    ASSERT_EQUAL("H2O", x.entity);
    ASSERT_EQUAL(1, x.value);
    ASSERT_EQUAL("kg", x.units);
}

auto test_MixtureParser() -> void
{
    const std::string s1 = "Mixture: 1 kg H2O; 1 mmol NaCl";
    const std::string s2 = "Mixture: \n - 1 kg H2O \n - 1 mmol NaCl";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);

    Mixture x1; node1 >> x1;
    Mixture x2; node2 >> x2;

    ASSERT_EQUAL("H2O", x1[0].entity);
    ASSERT_EQUAL("H2O", x2[0].entity);
    ASSERT_EQUAL(1, x1[0].value);
    ASSERT_EQUAL(1, x2[0].value);
    ASSERT_EQUAL("kg", x1[0].units);
    ASSERT_EQUAL("kg", x2[0].units);

    ASSERT_EQUAL("NaCl", x1[1].entity);
    ASSERT_EQUAL("NaCl", x2[1].entity);
    ASSERT_EQUAL(1, x1[1].value);
    ASSERT_EQUAL(1, x2[1].value);
    ASSERT_EQUAL("mmol", x1[1].units);
    ASSERT_EQUAL("mmol", x2[1].units);
}

auto test_EquilibriumConstraint_pH_Parser() -> void
{
    const std::string s1 = "pH: 3.0";
    const std::string s2 = "pH: 4.0 CO2";
    const std::string s3 = "pH: 5.0 HCl NaOH";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);

    EquilibriumConstraint::pH x1; node1 >> x1;
    EquilibriumConstraint::pH x2; node2 >> x2;
    EquilibriumConstraint::pH x3; node3 >> x3;

    ASSERT_EQUAL(3.0, x1.value);
    ASSERT_EQUAL(4.0, x2.value);
    ASSERT_EQUAL(5.0, x3.value);

    ASSERT_EQUAL("", x1.titrant1);
    ASSERT_EQUAL("", x1.titrant2);

    ASSERT_EQUAL("CO2", x2.titrant1);
    ASSERT_EQUAL("", x2.titrant2);

    ASSERT_EQUAL("HCl", x3.titrant1);
    ASSERT_EQUAL("NaOH", x3.titrant2);
}

auto test_EquilibriumConstraint_SpeciesAmount_Parser() -> void
{
    const std::string s1 = "SpeciesAmount: Calcite 100.0 g";
    const std::string s2 = "SpeciesAmount: CO2(g) 4.0 mol CO2";
    const std::string s3 = "SpeciesAmount: CO2(g) 4.0";
    const std::string s4 = "SpeciesAmount: ";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    EquilibriumConstraint::SpeciesAmount x1; node1 >> x1;
    EquilibriumConstraint::SpeciesAmount x2; node2 >> x2;
    EquilibriumConstraint::SpeciesAmount x3; ASSERT_THROWS(node3 >> x3, std::runtime_error);
    EquilibriumConstraint::SpeciesAmount x4; ASSERT_THROWS(node4 >> x4, std::runtime_error);

    ASSERT_EQUAL("Calcite", x1.entity);
    ASSERT_EQUAL("CO2(g)", x2.entity);

    ASSERT_EQUAL(100.0, x1.value);
    ASSERT_EQUAL(4.0, x2.value);

    ASSERT_EQUAL("g", x1.units);
    ASSERT_EQUAL("mol", x2.units);

    ASSERT_EQUAL("", x1.titrant1);
    ASSERT_EQUAL("", x1.titrant2);

    ASSERT_EQUAL("CO2", x2.titrant1);
    ASSERT_EQUAL("", x2.titrant2);
}

auto test_EquilibriumConstraint_SpeciesActivity_Parser() -> void
{
    const std::string s1 = "SpeciesActivity: O2(g) 40.0";
    const std::string s2 = "SpeciesActivity: H+ 1e-5 HCl NaOH";
    const std::string s3 = "SpeciesActivity: H+";
    const std::string s4 = "SpeciesActivity: 1.0";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    EquilibriumConstraint::SpeciesActivity x1; node1 >> x1;
    EquilibriumConstraint::SpeciesActivity x2; node2 >> x2;
    EquilibriumConstraint::SpeciesActivity x3; ASSERT_THROWS(node3 >> x3, std::runtime_error);
    EquilibriumConstraint::SpeciesActivity x4; ASSERT_THROWS(node4 >> x4, std::runtime_error);

    ASSERT_EQUAL("O2(g)", x1.entity);
    ASSERT_EQUAL("H+", x2.entity);

    ASSERT_EQUAL(40.0, x1.value);
    ASSERT_EQUAL(1e-5, x2.value);

    ASSERT_EQUAL("", x1.titrant1);
    ASSERT_EQUAL("", x1.titrant2);

    ASSERT_EQUAL("HCl", x2.titrant1);
    ASSERT_EQUAL("NaOH", x2.titrant2);
}

auto test_EquilibriumConstraint_PhaseAmount_Parser() -> void
{
    const std::string s1 = "PhaseAmount: Calcite 1.0 kg";
    const std::string s2 = "PhaseAmount: Aqueous 100 mol (1:kg:H2O)(1:mol:CaCl2)";
    const std::string s3 = "PhaseAmount: Gaseous 1.0";
    const std::string s4 = "PhaseAmount: ";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    EquilibriumConstraint::PhaseAmount x1; node1 >> x1;
    EquilibriumConstraint::PhaseAmount x2; node2 >> x2;
    EquilibriumConstraint::PhaseAmount x3; ASSERT_THROWS(node3 >> x3, std::runtime_error);
    EquilibriumConstraint::PhaseAmount x4; ASSERT_THROWS(node4 >> x4, std::runtime_error);

    ASSERT_EQUAL("Calcite", x1.entity);
    ASSERT_EQUAL("Aqueous", x2.entity);

    ASSERT_EQUAL(1.0, x1.value);
    ASSERT_EQUAL(100, x2.value);

    ASSERT_EQUAL("kg", x1.units);
    ASSERT_EQUAL("mol", x2.units);

    ASSERT_EQUAL("", x1.titrant1);
    ASSERT_EQUAL("", x1.titrant2);

    ASSERT_EQUAL("(1:kg:H2O)(1:mol:CaCl2)", x2.titrant1);
    ASSERT_EQUAL("", x2.titrant2);
}

auto test_EquilibriumConstraint_PhaseVolume_Parser() -> void
{
    const std::string s1 = "PhaseVolume: Calcite 1.0 m3";
    const std::string s2 = "PhaseVolume: Aqueous 100 cm3 (1:kg:H2O)(1:mol:CaCl2)";
    const std::string s3 = "PhaseVolume: Gaseous 1.0";
    const std::string s4 = "PhaseVolume: ";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    EquilibriumConstraint::PhaseVolume x1; node1 >> x1;
    EquilibriumConstraint::PhaseVolume x2; node2 >> x2;
    EquilibriumConstraint::PhaseVolume x3; ASSERT_THROWS(node3 >> x3, std::runtime_error);
    EquilibriumConstraint::PhaseVolume x4; ASSERT_THROWS(node4 >> x4, std::runtime_error);

    ASSERT_EQUAL("Calcite", x1.entity);
    ASSERT_EQUAL("Aqueous", x2.entity);

    ASSERT_EQUAL(1.0, x1.value);
    ASSERT_EQUAL(100, x2.value);

    ASSERT_EQUAL("m3", x1.units);
    ASSERT_EQUAL("cm3", x2.units);

    ASSERT_EQUAL("", x1.titrant1);
    ASSERT_EQUAL("", x1.titrant2);

    ASSERT_EQUAL("(1:kg:H2O)(1:mol:CaCl2)", x2.titrant1);
    ASSERT_EQUAL("", x2.titrant2);
}

auto test_EquilibriumConstraintEquilibriumParser() -> void
{
    std::string s1 = R"xyz(
Equilibrium: 
   Temperature: 30 celsius
   Pressure: 10 bar
   Mixture: 1 kg H2O; 1 mmol NaCl
)xyz";

    std::string s2 = R"xyz(
Equilibrium StateIC: 
   Temperature: 400 kelvin
   Pressure: 100 bar
   Mixture: 
     1 kg H2O
     1 mmol NaCl
   pH: 5.0 HCl
   SpeciesAmount: H2O(g) 1 ug
   Amount: Calcite 100 g
   Activity: O2(g) 0.20
   SpeciesActivity: CO2(g) 30 CO2
   PhaseVolume: Aqueous 1 m3 (1:kg:H2O)(1:mmol:NaCl)
   PhaseAmount: Magnesite 1 m3
   InertSpecies: Dolomite 100 g
   InertSpecies: Quartz 5 mg
   InertSpecies: Siderite 1 mg
   InertPhases: Gaseous Halite
)xyz";

    s1 = preprocess(s1);
    s2 = preprocess(s2);

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);

    Equilibrium x1; node1[0] >> x1;
    Equilibrium x2; node2[0] >> x2;

    ASSERT_EQUAL("", x1.stateid);

    ASSERT_EQUAL(30, x1.temperature.value);
    ASSERT_EQUAL("celsius", x1.temperature.units);

    ASSERT_EQUAL(10, x1.pressure.value);
    ASSERT_EQUAL("bar", x1.pressure.units);

    ASSERT_EQUAL(2, x1.mixture.size());
    ASSERT_EQUAL(1, x1.mixture[0].value);
    ASSERT_EQUAL("kg", x1.mixture[0].units);
    ASSERT_EQUAL("H2O", x1.mixture[0].entity);
    ASSERT_EQUAL(1, x1.mixture[1].value);
    ASSERT_EQUAL("mmol", x1.mixture[1].units);
    ASSERT_EQUAL("NaCl", x1.mixture[1].entity);


    ASSERT_EQUAL("StateIC", x2.stateid);

    ASSERT_EQUAL(400, x2.temperature.value);
    ASSERT_EQUAL("kelvin", x2.temperature.units);

    ASSERT_EQUAL(100, x2.pressure.value);
    ASSERT_EQUAL("bar", x2.pressure.units);

    ASSERT_EQUAL(2, x2.mixture.size());
    ASSERT_EQUAL(1, x2.mixture[0].value);
    ASSERT_EQUAL("kg", x2.mixture[0].units);
    ASSERT_EQUAL("H2O", x2.mixture[0].entity);
    ASSERT_EQUAL(1, x2.mixture[1].value);
    ASSERT_EQUAL("mmol", x2.mixture[1].units);
    ASSERT_EQUAL("NaCl", x2.mixture[1].entity);

    ASSERT_EQUAL(1, x2.pH.size());
    ASSERT_EQUAL(5.0, x2.pH[0].value);
    ASSERT_EQUAL("HCl", x2.pH[0].titrant1);

    ASSERT_EQUAL(2, x2.species_amounts.size());
    ASSERT_EQUAL("H2O(g)", x2.species_amounts[0].entity);
    ASSERT_EQUAL("ug", x2.species_amounts[0].units);
    ASSERT_EQUAL(1, x2.species_amounts[0].value);
    ASSERT_EQUAL("Calcite", x2.species_amounts[1].entity);
    ASSERT_EQUAL("g", x2.species_amounts[1].units);
    ASSERT_EQUAL(100, x2.species_amounts[1].value);

    ASSERT_EQUAL(2, x2.species_activities.size());
    ASSERT_EQUAL("O2(g)", x2.species_activities[0].entity);
    ASSERT_EQUAL(0.2, x2.species_activities[0].value);
    ASSERT_EQUAL("CO2(g)", x2.species_activities[1].entity);
    ASSERT_EQUAL(30, x2.species_activities[1].value);
    ASSERT_EQUAL("CO2", x2.species_activities[1].titrant1);

    ASSERT_EQUAL(1, x2.phase_volumes.size());
    ASSERT_EQUAL("Aqueous", x2.phase_volumes[0].entity);
    ASSERT_EQUAL(1, x2.phase_volumes[0].value);
    ASSERT_EQUAL("m3", x2.phase_volumes[0].units);
    ASSERT_EQUAL("(1:kg:H2O)(1:mmol:NaCl)", x2.phase_volumes[0].titrant1);

    ASSERT_EQUAL(1, x2.phase_amounts.size());
    ASSERT_EQUAL("Magnesite", x2.phase_amounts[0].entity);
    ASSERT_EQUAL(1, x2.phase_amounts[0].value);
    ASSERT_EQUAL("m3", x2.phase_amounts[0].units);

    ASSERT_EQUAL(3, x2.inert_species.size());
    ASSERT_EQUAL("Dolomite", x2.inert_species[0].entity);
    ASSERT_EQUAL(100, x2.inert_species[0].value);
    ASSERT_EQUAL("g", x2.inert_species[0].units);
    ASSERT_EQUAL("Quartz", x2.inert_species[1].entity);
    ASSERT_EQUAL(5, x2.inert_species[1].value);
    ASSERT_EQUAL("mg", x2.inert_species[1].units);
    ASSERT_EQUAL("Siderite", x2.inert_species[2].entity);
    ASSERT_EQUAL(1, x2.inert_species[2].value);
    ASSERT_EQUAL("mg", x2.inert_species[2].units);

    ASSERT_EQUAL(2, x2.inert_phases.size());
    ASSERT_EQUAL("Gaseous", x2.inert_phases[0]);
    ASSERT_EQUAL("Halite", x2.inert_phases[1]);
}

int main()
{
    cute::suite s;
    s += CUTE(test_ValueUnitsParser);
    s += CUTE(test_EntityValueUnitsParser);
    s += CUTE(test_MixtureCompoundParser);
    s += CUTE(test_MixtureParser);
    s += CUTE(test_EquilibriumConstraint_pH_Parser);
    s += CUTE(test_EquilibriumConstraint_SpeciesAmount_Parser);
    s += CUTE(test_EquilibriumConstraint_SpeciesActivity_Parser);
    s += CUTE(test_EquilibriumConstraint_PhaseAmount_Parser);
    s += CUTE(test_EquilibriumConstraint_PhaseVolume_Parser);
    s += CUTE(test_EquilibriumConstraintEquilibriumParser);

    cute::ide_listener<> lis;
    cute::makeRunner(lis)(s, "TestParseUtils");
}
