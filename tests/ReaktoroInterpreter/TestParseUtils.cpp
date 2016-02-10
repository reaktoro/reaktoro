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
using namespace iReaktoro;

// cute includes
#include <cute/cute.h>
#include <cute/cute_runner.h>
#include <cute/ide_listener.h>

auto test_ValueUnitsParser() -> void
{
    const std::string s = "Temperature: 25 celsius";
    YAML::Node node = YAML::Load(s);
    ValueUnits x; node["Temperature"] >> x;

    ASSERT_EQUAL(25, x.value);
    ASSERT_EQUAL("celsius", x.units);
}

auto test_EntityValueUnitsParser() -> void
{
    const std::string s = "SpeciesAmount: Calcite 100 g";
    YAML::Node node = YAML::Load(s);
    EntityValueUnits x; node["SpeciesAmount"] >> x;

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

    Mixture x1; node1["Mixture"] >> x1;
    Mixture x2; node2["Mixture"] >> x2;

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

    EquilibriumConstraint::pH x1; node1["pH"] >> x1;
    EquilibriumConstraint::pH x2; node2["pH"] >> x2;
    EquilibriumConstraint::pH x3; node3["pH"] >> x3;

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

    EquilibriumConstraint::SpeciesAmount x1; node1["SpeciesAmount"] >> x1;
    EquilibriumConstraint::SpeciesAmount x2; node2["SpeciesAmount"] >> x2;
    EquilibriumConstraint::SpeciesAmount x3; ASSERT_THROWS(node3["SpeciesAmount"] >> x3, std::runtime_error);
    EquilibriumConstraint::SpeciesAmount x4; ASSERT_THROWS(node4["SpeciesAmount"] >> x4, std::runtime_error);

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

    EquilibriumConstraint::SpeciesAmount x1; node1["SpeciesActivity"] >> x1;
    EquilibriumConstraint::SpeciesAmount x2; node2["SpeciesActivity"] >> x2;
    EquilibriumConstraint::SpeciesAmount x3; ASSERT_THROWS(node3["SpeciesActivity"] >> x3, std::runtime_error);
    EquilibriumConstraint::SpeciesAmount x4; ASSERT_THROWS(node4["SpeciesActivity"] >> x4, std::runtime_error);

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

    EquilibriumConstraint::SpeciesAmount x1; node1["PhaseAmount"] >> x1;
    EquilibriumConstraint::SpeciesAmount x2; node2["PhaseAmount"] >> x2;
    EquilibriumConstraint::SpeciesAmount x3; ASSERT_THROWS(node3["PhaseAmount"] >> x3, std::runtime_error);
    EquilibriumConstraint::SpeciesAmount x4; ASSERT_THROWS(node4["PhaseAmount"] >> x4, std::runtime_error);

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
    const std::string s1 = "PhaseAmount: Calcite 1.0 m3";
    const std::string s2 = "PhaseAmount: Aqueous 100 cm3 (1:kg:H2O)(1:mol:CaCl2)";
    const std::string s3 = "PhaseAmount: Gaseous 1.0";
    const std::string s4 = "PhaseAmount: ";

    YAML::Node node1 = YAML::Load(s1);
    YAML::Node node2 = YAML::Load(s2);
    YAML::Node node3 = YAML::Load(s3);
    YAML::Node node4 = YAML::Load(s4);

    EquilibriumConstraint::SpeciesAmount x1; node1["PhaseAmount"] >> x1;
    EquilibriumConstraint::SpeciesAmount x2; node2["PhaseAmount"] >> x2;
    EquilibriumConstraint::SpeciesAmount x3; ASSERT_THROWS(node3["PhaseAmount"] >> x3, std::runtime_error);
    EquilibriumConstraint::SpeciesAmount x4; ASSERT_THROWS(node4["PhaseAmount"] >> x4, std::runtime_error);

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
    cute::ide_listener<> lis;
    cute::makeRunner(lis)(s, "TestParseUtils");
}
