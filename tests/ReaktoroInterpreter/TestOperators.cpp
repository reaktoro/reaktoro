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

// Reaktoro includes
#include <Reaktoro/Interpreter/Operators.hpp>
using namespace Reaktoro::kwd;

// cute includes
#include <cute/cute.h>
#include <cute/cute_runner.h>
#include <cute/ide_listener.h>

auto test_ShiftOperatorValueUnits() -> void
{
    const std::string s = "Temperature: 25 celsius";
    YAML::Node node = YAML::Load(s);
    ValueUnits x; node >> x;

    ASSERT_EQUAL(25, x.value);
    ASSERT_EQUAL("celsius", x.units);
}

auto test_ShiftOperatorEntityValueUnits() -> void
{
    const std::string s = "SpeciesAmount: Calcite 100 g";
    YAML::Node node = YAML::Load(s);
    EntityValueUnits x; node >> x;

    ASSERT_EQUAL("Calcite", x.entity);
    ASSERT_EQUAL(100, x.value);
    ASSERT_EQUAL("g", x.units);
}

auto test_ShiftOperatorValueUnitsEntity() -> void
{
    const std::string s = "1 kg H2O";
    YAML::Node node = YAML::Load(s);
    ValueUnitsEntity x; node >> x;

    ASSERT_EQUAL("H2O", x.entity);
    ASSERT_EQUAL(1, x.value);
    ASSERT_EQUAL("kg", x.units);
}

auto test_ShiftOperatorRecipe() -> void
{
    const std::string s1 = "Recipe: 1 kg H2O; 1 mmol NaCl";
    const std::string s2 = "Recipe: \n - 1 kg H2O \n - 1 mmol NaCl";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);

    Recipe x1; node1 >> x1;
    Recipe x2; node2 >> x2;

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

auto test_ShiftOperatorEquilibriumConstraint_pH() -> void
{
    const std::string s1 = "pH: 3.0";
    const std::string s2 = "pH: 4.0 CO2";
    const std::string s3 = "pH: 5.0 HCl NaOH";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);

    pH x1; node1 >> x1;
    pH x2; node2 >> x2;
    pH x3; node3 >> x3;

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

auto test_ShiftOperatorEquilibriumConstraint_SpeciesAmount() -> void
{
    const std::string s1 = "SpeciesAmount: Calcite 100.0 g";
    const std::string s2 = "SpeciesAmount: CO2(g) 4.0 mol CO2";
    const std::string s3 = "SpeciesAmount: CO2(g) 4.0";
    const std::string s4 = "SpeciesAmount: ";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    SpeciesAmount x1; node1 >> x1;
    SpeciesAmount x2; node2 >> x2;
    SpeciesAmount x3; ASSERT_THROWS(node3 >> x3, std::runtime_error);
    SpeciesAmount x4; ASSERT_THROWS(node4 >> x4, std::runtime_error);

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

auto test_ShiftOperatorEquilibriumConstraint_SpeciesActivity() -> void
{
    const std::string s1 = "SpeciesActivity: O2(g) 40.0";
    const std::string s2 = "SpeciesActivity: H+ 1e-5 HCl NaOH";
    const std::string s3 = "SpeciesActivity: H+";
    const std::string s4 = "SpeciesActivity: 1.0";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    SpeciesActivity x1; node1 >> x1;
    SpeciesActivity x2; node2 >> x2;
    SpeciesActivity x3; ASSERT_THROWS(node3 >> x3, std::runtime_error);
    SpeciesActivity x4; ASSERT_THROWS(node4 >> x4, std::runtime_error);

    ASSERT_EQUAL("O2(g)", x1.entity);
    ASSERT_EQUAL("H+", x2.entity);

    ASSERT_EQUAL(40.0, x1.value);
    ASSERT_EQUAL(1e-5, x2.value);

    ASSERT_EQUAL("", x1.titrant1);
    ASSERT_EQUAL("", x1.titrant2);

    ASSERT_EQUAL("HCl", x2.titrant1);
    ASSERT_EQUAL("NaOH", x2.titrant2);
}

auto test_ShiftOperatorEquilibriumConstraint_PhaseAmount() -> void
{
    const std::string s1 = "PhaseAmount: Calcite 1.0 kg";
    const std::string s2 = "PhaseAmount: Aqueous 100 mol (1:kg:H2O)(1:mol:CaCl2)";
    const std::string s3 = "PhaseAmount: Gaseous 1.0";
    const std::string s4 = "PhaseAmount: ";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    PhaseAmount x1; node1 >> x1;
    PhaseAmount x2; node2 >> x2;
    PhaseAmount x3; ASSERT_THROWS(node3 >> x3, std::runtime_error);
    PhaseAmount x4; ASSERT_THROWS(node4 >> x4, std::runtime_error);

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

auto test_ShiftOperatorEquilibriumConstraint_PhaseVolume() -> void
{
    const std::string s1 = "PhaseVolume: Calcite 1.0 m3";
    const std::string s2 = "PhaseVolume: Aqueous 100 cm3 (1:kg:H2O)(1:mol:CaCl2)";
    const std::string s3 = "PhaseVolume: Gaseous 1.0";
    const std::string s4 = "PhaseVolume: ";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    PhaseVolume x1; node1 >> x1;
    PhaseVolume x2; node2 >> x2;
    PhaseVolume x3; ASSERT_THROWS(node3 >> x3, std::runtime_error);
    PhaseVolume x4; ASSERT_THROWS(node4 >> x4, std::runtime_error);

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

auto test_ShiftOperatorPlot() -> void
{
    std::string s = R"xyz(
- Plot:
   - Name: Calcite
   - x: pH
   - y: molality species=Ca++ units=mmolal
   - xlabel: pH
   - ylabel: Concentration [mmolal]
   - ytitles: Ca++
   - Key: right center
)xyz";

    YAML::Node node = YAML::Load(s);

    Plot x; node[0] >> x;

    ASSERT_EQUAL("Calcite", x.name);
    ASSERT_EQUAL("pH", x.x);
    ASSERT_EQUAL("molality species=Ca++ units=mmolal", x.y);
    ASSERT_EQUAL("pH", x.xlabel);
    ASSERT_EQUAL("Concentration [mmolal]", x.ylabel);
    ASSERT_EQUAL("Ca++", x.ytitles);
    ASSERT_EQUAL("right center", x.key);
}

auto test_ShiftOperatorEquilibrium() -> void
{
    std::string s1 = R"xyz(
- Equilibrium:
   - Temperature: 30 celsius
   - Pressure: 10 bar
   - Recipe: 1 kg H2O; 1 mmol NaCl
)xyz";

    std::string s2 = R"xyz(
- Equilibrium StateIC:
   - Temperature: 400 kelvin
   - Pressure: 100 bar
   - Recipe:
     - 1 kg H2O
     - 1 mmol NaCl
   - pH: 5.0 HCl
   - SpeciesAmount: H2O(g) 1 ug
   - Amount: Calcite 100 g
   - Activity: O2(g) 0.20
   - SpeciesActivity: CO2(g) 30 CO2
   - PhaseVolume: Aqueous 1 m3 (1:kg:H2O)(1:mmol:NaCl)
   - PhaseAmount: Magnesite 1 m3
   - InertSpecies: Dolomite 100 g
   - InertSpecies: Quartz 5 mg
   - InertSpecies: Siderite 1 mg
   - InertPhases: Gaseous Halite
)xyz";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);

    EquilibriumProblem x1; node1[0] >> x1;
    EquilibriumProblem x2; node2[0] >> x2;

    ASSERT_EQUAL("", x1.stateid);

    ASSERT_EQUAL(30, x1.temperature.value);
    ASSERT_EQUAL("celsius", x1.temperature.units);

    ASSERT_EQUAL(10, x1.pressure.value);
    ASSERT_EQUAL("bar", x1.pressure.units);

    ASSERT_EQUAL(2, x1.recipe.size());
    ASSERT_EQUAL(1, x1.recipe[0].value);
    ASSERT_EQUAL("kg", x1.recipe[0].units);
    ASSERT_EQUAL("H2O", x1.recipe[0].entity);
    ASSERT_EQUAL(1, x1.recipe[1].value);
    ASSERT_EQUAL("mmol", x1.recipe[1].units);
    ASSERT_EQUAL("NaCl", x1.recipe[1].entity);


    ASSERT_EQUAL("StateIC", x2.stateid);

    ASSERT_EQUAL(400, x2.temperature.value);
    ASSERT_EQUAL("kelvin", x2.temperature.units);

    ASSERT_EQUAL(100, x2.pressure.value);
    ASSERT_EQUAL("bar", x2.pressure.units);

    ASSERT_EQUAL(2, x2.recipe.size());
    ASSERT_EQUAL(1, x2.recipe[0].value);
    ASSERT_EQUAL("kg", x2.recipe[0].units);
    ASSERT_EQUAL("H2O", x2.recipe[0].entity);
    ASSERT_EQUAL(1, x2.recipe[1].value);
    ASSERT_EQUAL("mmol", x2.recipe[1].units);
    ASSERT_EQUAL("NaCl", x2.recipe[1].entity);

    ASSERT_EQUAL(1, x2.ph.size());
    ASSERT_EQUAL(5.0, x2.ph[0].value);
    ASSERT_EQUAL("HCl", x2.ph[0].titrant1);

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

auto test_ShiftOperatorEquilibriumPath() -> void
{
    std::string s = R"xyz(
- EquilibriumPath:
   - InitialState: State1
   - FinalState: State2
   - Plot:
       - Name: Calcite
       - x: t
       - y: mass species=Calcite units=g
       - xlabel: t [hour]
       - ylabel: Concentration [mmolal]
       - ytitles: Calcite
       - Key: right center
)xyz";

    YAML::Node node = YAML::Load(s);

    EquilibriumPath x; node[0] >> x;

    ASSERT_EQUAL("State1", x.initial_state);
    ASSERT_EQUAL("State2", x.final_state);
    ASSERT_EQUAL(1, x.plots.size());
    ASSERT_EQUAL("Calcite", x.plots[0].name);
}

auto test_ShiftOperatorKineticPath() -> void
{
    std::string s = R"xyz(
- KineticPath StateFinal:
   - InitialCondition: StateIC
   - KineticSpecies: Calcite Dolomite
   - Duration: 24 hours
   - Plot:
       - Name: MolalityCa
       - x: t
       - y: molality element=Ca units=mmolal
       - xlabel: t [hour]
       - ylabel: Concentration [mmolal]
       - ytitles: Ca
       - Key: right center
   - Plot:
       - Name: Calcite
       - x: t
       - y: mass species=Calcite units=g
       - xlabel: t [hour]
       - ylabel: Concentration [mmolal]
       - ytitles: Calcite
       - Key: right center
)xyz";

    YAML::Node node = YAML::Load(s);

    KineticPath x; node[0] >> x;

    ASSERT_EQUAL("StateFinal", x.stateid);
    ASSERT_EQUAL("StateIC", x.initial_condition);
    ASSERT_EQUAL(2, x.kinetic_species.size());
    ASSERT_EQUAL("Calcite", x.kinetic_species[0]);
    ASSERT_EQUAL("Dolomite", x.kinetic_species[1]);
    ASSERT_EQUAL(2, x.duration.value);
    ASSERT_EQUAL("hours", x.duration.units);
    ASSERT_EQUAL(2, x.plots.size());
    ASSERT_EQUAL("MolalityCa", x.plots[0].name);
    ASSERT_EQUAL("Calcite", x.plots[1].name);
}

auto test_ShiftOperatorMineralReaction() -> void
{
    std::string s = R"xyz(
- MineralReaction Calcite:
    - Equation: Calcite = Ca++ + CO3--
    - Mechanism: logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol
    - Mechanism: logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0
    - SpecificSurfaceArea: 10 cm2/g
)xyz";

    YAML::Node node = YAML::Load(s);

    MineralReaction x; node[0] >> x;

    ASSERT_EQUAL("Calcite", x.mineral);
    ASSERT_EQUAL("Calcite = Ca++ + CO3--", x.equation);
    ASSERT_EQUAL(2, x.mechanisms.size());
    ASSERT_EQUAL("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol", x.mechanisms[0]);
    ASSERT_EQUAL("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0", x.mechanisms[1]);
    ASSERT_EQUAL(10, x.ssa.value);
    ASSERT_EQUAL("cm2/g", x.ssa.units);
}

auto test_ShiftOperatorConcentrations() -> void
{
    const std::string s1 = "Concentrations: Na 200; Cl 200";
    const std::string s2 = "Concentrations:\n  Na 200\n  Cl 200";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);

    Concentrations x1; node1 >> x1;
    Concentrations x2; node2 >> x2;

    ASSERT_EQUAL("Na", x1[0].entity);
    ASSERT_EQUAL("Na", x2[0].entity);
    ASSERT_EQUAL(200, x1[0].value);
    ASSERT_EQUAL(200, x2[0].value);

    ASSERT_EQUAL("Cl", x1[1].entity);
    ASSERT_EQUAL("Cl", x2[1].entity);
    ASSERT_EQUAL(200, x1[1].value);
    ASSERT_EQUAL(200, x2[1].value);
}

auto test_ShiftOperatorSpeciationProblem() -> void
{
    std::string s1 = R"xyz(
- Speciation:
   - Temperature: 30 celsius
   - Pressure: 10 bar
   - Concentrations: Na 200; Cl 200
)xyz";

    std::string s2 = R"xyz(
- Speciation StateSpec:
   - Temperature: 400 kelvin
   - Pressure: 100 bar
   - Concentrations:
     - Na  200
     - Cl  200
   - pH: 5.0 HCl
   - Fugacity: O2(g) 4 bar
   - InertSpecies: Dolomite 100 g
   - InertSpecies: Quartz 5 mg
   - InertSpecies: Siderite 1 mg
   - InertPhases: Gaseous Halite
)xyz";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);

    SpeciationProblem x1; node1[0] >> x1;
    SpeciationProblem x2; node2[0] >> x2;

    ASSERT_EQUAL("", x1.stateid);

    ASSERT_EQUAL(30, x1.temperature.value);
    ASSERT_EQUAL("celsius", x1.temperature.units);

    ASSERT_EQUAL(10, x1.pressure.value);
    ASSERT_EQUAL("bar", x1.pressure.units);

    ASSERT_EQUAL(2, x1.concentrations.size());
    ASSERT_EQUAL("Na", x1.concentrations[0].entity);
    ASSERT_EQUAL(200, x1.concentrations[0].value);
    ASSERT_EQUAL("Cl", x1.concentrations[1].entity);
    ASSERT_EQUAL(200, x1.concentrations[1].value);


    ASSERT_EQUAL("StateSpec", x2.stateid);

    ASSERT_EQUAL(400, x2.temperature.value);
    ASSERT_EQUAL("kelvin", x2.temperature.units);

    ASSERT_EQUAL(100, x2.pressure.value);
    ASSERT_EQUAL("bar", x2.pressure.units);

    ASSERT_EQUAL(2, x2.concentrations.size());
    ASSERT_EQUAL("Na", x2.concentrations[0].entity);
    ASSERT_EQUAL(200, x2.concentrations[0].value);
    ASSERT_EQUAL("Cl", x2.concentrations[1].entity);
    ASSERT_EQUAL(200, x2.concentrations[1].value);

    ASSERT_EQUAL(1, x2.ph.size());
    ASSERT_EQUAL(5.0, x2.ph[0].value);
    ASSERT_EQUAL("HCl", x2.ph[0].titrant1);

    ASSERT_EQUAL(1, x2.species_fugacities.size());
    ASSERT_EQUAL("O2(g)", x2.species_fugacities[0].entity);
    ASSERT_EQUAL(4, x2.species_fugacities[0].value);
    ASSERT_EQUAL("bar", x2.species_fugacities[0].units);

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
auto test_ShiftOperatorPhreeqcKeyword() -> void
{
    std::string s1 = R"xyz(
- PHREEQC StatePHREEQC:
   - Database: pitzer.dat
   - Input: input.dat
)xyz";

    std::string s2 = R"xyz(
- PHREEQC:
   - Database: phreeqc.dat
   - Input:
      - SOLUTION
      - EQUILIBRIUM_PHASES
)xyz";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);

    PhreeqcKeyword x1; node1[0] >> x1;
    PhreeqcKeyword x2; node2[0] >> x2;

    ASSERT_EQUAL("StatePHREEQC", x1.stateid);
    ASSERT_EQUAL("StatePhreeqc", x2.stateid);

    ASSERT_EQUAL("pitzer.dat", x1.database);
    ASSERT_EQUAL("phreeqc.dat", x2.database);

    ASSERT_EQUAL("input.dat", x1.input);
    ASSERT_EQUAL("SOLUTION\nEQUILIBRIUM_PHASES", x2.input);
}

int main()
{
    cute::suite s;
    s += CUTE(test_ShiftOperatorValueUnits);
    s += CUTE(test_ShiftOperatorEntityValueUnits);
    s += CUTE(test_ShiftOperatorValueUnitsEntity);
    s += CUTE(test_ShiftOperatorRecipe);
    s += CUTE(test_ShiftOperatorEquilibriumConstraint_pH);
    s += CUTE(test_ShiftOperatorEquilibriumConstraint_SpeciesAmount);
    s += CUTE(test_ShiftOperatorEquilibriumConstraint_SpeciesActivity);
    s += CUTE(test_ShiftOperatorEquilibriumConstraint_PhaseAmount);
    s += CUTE(test_ShiftOperatorEquilibriumConstraint_PhaseVolume);
    s += CUTE(test_ShiftOperatorPlot);
    s += CUTE(test_ShiftOperatorEquilibrium);
    s += CUTE(test_ShiftOperatorEquilibriumPath);
    s += CUTE(test_ShiftOperatorKineticPath);
    s += CUTE(test_ShiftOperatorConcentrations);
    s += CUTE(test_ShiftOperatorSpeciationProblem);
    s += CUTE(test_ShiftOperatorPhreeqcKeyword);

    cute::ide_listener<> lis;
    cute::makeRunner(lis)(s, "TestParseUtils");
}
