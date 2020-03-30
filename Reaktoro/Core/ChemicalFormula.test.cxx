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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalFormula.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ChemicalFormula class", "[ChemicalFormula]")
{
    ChemicalFormula formula;

    formula = ChemicalFormula("H2O");
    REQUIRE(formula.str() == "H2O");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 2);
    REQUIRE(formula.coefficient("H") == 2);
    REQUIRE(formula.coefficient("O") == 1);
    REQUIRE(formula.equivalent("H2O(aq)"));
    REQUIRE(formula.equivalent("H2O(l)"));
    REQUIRE(formula.equivalent("HHO"));
    REQUIRE(formula.equivalent("HOH"));
    REQUIRE(formula.equivalent("OH2"));

    formula = ChemicalFormula("CaCO3");
    REQUIRE(formula.str() == "CaCO3");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 3);
    REQUIRE(formula.coefficient("C") == 1);
    REQUIRE(formula.coefficient("Ca") == 1);
    REQUIRE(formula.coefficient("O") == 3);
    REQUIRE(formula.equivalent("CaCOOO"));
    REQUIRE(formula.equivalent("CaOOOC"));
    REQUIRE(formula.equivalent("Ca(CO3)"));
    REQUIRE(formula.equivalent("Ca(CO3)(aq)"));
    REQUIRE(formula.equivalent("Ca(CO3)(s)"));

    formula = ChemicalFormula("HCO3-");
    REQUIRE(formula.str() == "HCO3-");
    REQUIRE(formula.charge() == -1);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 4);
    REQUIRE(formula.coefficient("C") == 1);
    REQUIRE(formula.coefficient("H") == 1);
    REQUIRE(formula.coefficient("O") == 3);
    REQUIRE(formula.equivalent("HCO3-(aq)"));
    REQUIRE(formula.equivalent("HCO3(-)(aq)"));
    REQUIRE(formula.equivalent("HCO3(-)"));
    REQUIRE(formula.equivalent("HCOOO-"));

    formula = ChemicalFormula("H+");
    REQUIRE(formula.str() == "H+");
    REQUIRE(formula.charge() == 1);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 2);
    REQUIRE(formula.coefficient("H") == 1);
    REQUIRE(formula.equivalent("H+(aq)"));
    REQUIRE(formula.equivalent("H(+)"));

    formula = ChemicalFormula("Na+");
    REQUIRE(formula.str() == "Na+");
    REQUIRE(formula.charge() == 1);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 2);
    REQUIRE(formula.coefficient("Na") == 1);
    REQUIRE(formula.equivalent("Na+(aq)"));
    REQUIRE(formula.equivalent("Na+(pl)"));
    REQUIRE(formula.equivalent("Na(+)(aq)"));
    REQUIRE(formula.equivalent("Na(+)"));

    formula = ChemicalFormula("Cl-");
    REQUIRE(formula.str() == "Cl-");
    REQUIRE(formula.charge() == -1);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 2);
    REQUIRE(formula.coefficient("Cl") == 1);
    REQUIRE(formula.equivalent("Cl(-)(aq)"));
    REQUIRE(formula.equivalent("Cl-(aq)"));
    REQUIRE(formula.equivalent("Cl-(pl)"));

    formula = ChemicalFormula("CO3--");
    REQUIRE(formula.str() == "CO3--");
    REQUIRE(formula.charge() == -2);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 3);
    REQUIRE(formula.coefficient("C") == 1);
    REQUIRE(formula.coefficient("O") == 3);
    REQUIRE(formula.equivalent("CO3-2(aq)"));
    REQUIRE(formula.equivalent("CO3-2"));
    REQUIRE(formula.equivalent("CO3(2-)(aq)"));

    formula = ChemicalFormula("Fe+++");
    REQUIRE(formula.str() == "Fe+++");
    REQUIRE(formula.charge() == 3);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 2);
    REQUIRE(formula.coefficient("Fe") == 1);
    REQUIRE(formula.equivalent("Fe+++(aq)"));
    REQUIRE(formula.equivalent("Fe+3(aq)"));
    REQUIRE(formula.equivalent("Fe(3+)(aq)"));
    REQUIRE(formula.equivalent("Fe+3"));
    REQUIRE(formula.equivalent("Fe(3+)"));

    formula = ChemicalFormula("(CaMg)(CO3)2");
    REQUIRE(formula.str() == "(CaMg)(CO3)2");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 4);
    REQUIRE(formula.coefficient("C") == 2);
    REQUIRE(formula.coefficient("Ca") == 1);
    REQUIRE(formula.coefficient("Mg") == 1);
    REQUIRE(formula.coefficient("O") == 6);
    REQUIRE(formula.equivalent("CaMg(CO3)2(s)"));
    REQUIRE(formula.equivalent("CaMg(CO3)2"));
    REQUIRE(formula.equivalent("CaMgCO3CO3(s)"));
    REQUIRE(formula.equivalent("CaMgCO3CO3"));

    formula = ChemicalFormula("CH3COOH");
    REQUIRE(formula.str() == "CH3COOH");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 3);
    REQUIRE(formula.coefficient("C") == 2);
    REQUIRE(formula.coefficient("H") == 4);
    REQUIRE(formula.coefficient("O") == 2);
    REQUIRE(formula.equivalent("CH3COOH(aq)"));
    REQUIRE(formula.equivalent("CHHHCOOH(aq)"));
    REQUIRE(formula.equivalent("CHHHCOOH"));

    formula = ChemicalFormula("Al2.5Si0.5O4.75");
    REQUIRE(formula.str() == "Al2.5Si0.5O4.75");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 3);
    REQUIRE(formula.coefficient("Al") == 2.5);
    REQUIRE(formula.coefficient("Si") == 0.5);
    REQUIRE(formula.coefficient("O") == 4.75);
    REQUIRE(formula.equivalent("Al2.5Si0.5O4.75(s)"));
    REQUIRE(formula.equivalent("Si0.5Al2.5O4.75(s)"));

    formula = ChemicalFormula("Fe4Al18Si7.5O48H4");
    REQUIRE(formula.str() == "Fe4Al18Si7.5O48H4");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 5);
    REQUIRE(formula.coefficient("Fe") == 4);
    REQUIRE(formula.coefficient("Al") == 18);
    REQUIRE(formula.coefficient("Si") == 7.5);
    REQUIRE(formula.coefficient("O") == 48);
    REQUIRE(formula.coefficient("H") == 4);
    REQUIRE(formula.equivalent("Fe4Al18Si7.5O48H4(s)"));
    REQUIRE(formula.equivalent("Al18Fe4H4Si7.5O48(s)"));

    formula = ChemicalFormula("Mg4Al18Si7.5O48H4");
    REQUIRE(formula.str() == "Mg4Al18Si7.5O48H4");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 5);
    REQUIRE(formula.coefficient("Mg") == 4);
    REQUIRE(formula.coefficient("Al") == 18);
    REQUIRE(formula.coefficient("Si") == 7.5);
    REQUIRE(formula.coefficient("O") == 48);
    REQUIRE(formula.coefficient("H") == 4);
    REQUIRE(formula.equivalent("Mg4Al18Si7.5O48H4(s)"));

    formula = ChemicalFormula("Mn4Al18Si7.5O48H4");
    REQUIRE(formula.str() == "Mn4Al18Si7.5O48H4");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 5);
    REQUIRE(formula.coefficient("Mn") == 4);
    REQUIRE(formula.coefficient("Al") == 18);
    REQUIRE(formula.coefficient("Si") == 7.5);
    REQUIRE(formula.coefficient("O") == 48);
    REQUIRE(formula.coefficient("H") == 4);

    formula = ChemicalFormula("Ca0.5Al1Si2O6");
    REQUIRE(formula.str() == "Ca0.5Al1Si2O6");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 4);
    REQUIRE(formula.coefficient("Ca") == 0.5);
    REQUIRE(formula.coefficient("Al") == 1);
    REQUIRE(formula.coefficient("Si") == 2);
    REQUIRE(formula.coefficient("O") == 6);

    formula = ChemicalFormula("K0.5Fe5Al2Si8O30.5H12.5");
    REQUIRE(formula.str() == "K0.5Fe5Al2Si8O30.5H12.5");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 6);
    REQUIRE(formula.coefficient("K") == 0.5);
    REQUIRE(formula.coefficient("Fe") == 5);
    REQUIRE(formula.coefficient("Al") == 2);
    REQUIRE(formula.coefficient("Si") == 8);
    REQUIRE(formula.coefficient("O") == 30.5);
    REQUIRE(formula.coefficient("H") == 12.5);

    formula = ChemicalFormula("K0.5Mg5Al2Si8O30.5H12.5");
    REQUIRE(formula.str() == "K0.5Mg5Al2Si8O30.5H12.5");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 6);
    REQUIRE(formula.coefficient("K") == 0.5);
    REQUIRE(formula.coefficient("Mg") == 5);
    REQUIRE(formula.coefficient("Al") == 2);
    REQUIRE(formula.coefficient("Si") == 8);
    REQUIRE(formula.coefficient("O") == 30.5);
    REQUIRE(formula.coefficient("H") == 12.5);

    formula = ChemicalFormula("Mg3.5Al9Si1.5O20");
    REQUIRE(formula.str() == "Mg3.5Al9Si1.5O20");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 4);
    REQUIRE(formula.coefficient("Mg") == 3.5);
    REQUIRE(formula.coefficient("Al") == 9);
    REQUIRE(formula.coefficient("Si") == 1.5);
    REQUIRE(formula.coefficient("O") == 20);

    formula = ChemicalFormula("Fe3.5Al9Si1.5O20");
    REQUIRE(formula.str() == "Fe3.5Al9Si1.5O20");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 4);
    REQUIRE(formula.coefficient("Fe") == 3.5);
    REQUIRE(formula.coefficient("Al") == 9);
    REQUIRE(formula.coefficient("Si") == 1.5);
    REQUIRE(formula.coefficient("O") == 20);

    formula = ChemicalFormula("Fe0.875S1");
    REQUIRE(formula.str() == "Fe0.875S1");
    REQUIRE(formula.charge() == 0);
    REQUIRE(formula.molarMass() == Approx(0.0));
    REQUIRE(formula.symbols().size() == 2);
    REQUIRE(formula.coefficient("Fe") == 0.875);
    REQUIRE(formula.coefficient("S") == 1);

    REQUIRE(ChemicalFormula::equivalent("Ca++", "Ca+2"));
    REQUIRE(ChemicalFormula::equivalent("Ca++", "Ca(2+)"));
    REQUIRE(ChemicalFormula::equivalent("Ca++", "Ca(2+)(aq)"));

    REQUIRE(ChemicalFormula::equivalent("CO3--", "CO3-2"));
    REQUIRE(ChemicalFormula::equivalent("CO3--", "CO3(2-)"));
    REQUIRE(ChemicalFormula::equivalent("CO3--", "CO3(2-)(aq)"));

    REQUIRE(ChemicalFormula::equivalent("Fe+++", "Fe+3"));
    REQUIRE(ChemicalFormula::equivalent("Fe+++", "Fe(3+)"));
    REQUIRE(ChemicalFormula::equivalent("Fe+++", "Fe(3+)(aq)"));

    REQUIRE(ChemicalFormula::equivalent("H+", "H+1"));
    REQUIRE(ChemicalFormula::equivalent("H+", "H(+)"));
    REQUIRE(ChemicalFormula::equivalent("H+", "H(+)(aq)"));

    REQUIRE(ChemicalFormula::equivalent("OH-", "OH-1"));
    REQUIRE(ChemicalFormula::equivalent("OH-", "OH(-)"));
    REQUIRE(ChemicalFormula::equivalent("OH-", "OH(-)(aq)"));

    REQUIRE(ChemicalFormula::equivalent("CO2", "CO2(g)"));
    REQUIRE(ChemicalFormula::equivalent("CO2", "COO"));
}

TEST_CASE("Testing parseChemicalFormula method", "[ChemicalFormula]")
{
    std::unordered_map<std::string, double> formula;

    formula = parseChemicalFormula("H2O");
    REQUIRE(formula.size() == 2);
    REQUIRE(formula["H"] == 2);
    REQUIRE(formula["O"] == 1);

    formula = parseChemicalFormula("CaCO3");
    REQUIRE(formula.size() == 3);
    REQUIRE(formula["C"] == 1);
    REQUIRE(formula["Ca"] == 1);
    REQUIRE(formula["O"] == 3);

    formula = parseChemicalFormula("HCO3-");
    REQUIRE(formula.size() == 4);
    REQUIRE(formula["C"] == 1);
    REQUIRE(formula["H"] == 1);
    REQUIRE(formula["O"] == 3);
    REQUIRE(formula["Z"] == -1);

    formula = parseChemicalFormula("H+");
    REQUIRE(formula.size() == 2);
    REQUIRE(formula["H"] == 1);
    REQUIRE(formula["Z"] == 1);

    formula = parseChemicalFormula("Na+");
    REQUIRE(formula.size() == 2);
    REQUIRE(formula["Na"] == 1);
    REQUIRE(formula["Z"] == 1);

    formula = parseChemicalFormula("Cl-");
    REQUIRE(formula.size() == 2);
    REQUIRE(formula["Cl"] == 1);
    REQUIRE(formula["Z"] == -1);

    formula = parseChemicalFormula("CO3--");
    REQUIRE(formula.size() == 3);
    REQUIRE(formula["C"] == 1);
    REQUIRE(formula["O"] == 3);
    REQUIRE(formula["Z"] == -2);

    formula = parseChemicalFormula("CO3-2");
    REQUIRE(formula.size() == 3);
    REQUIRE(formula["C"] == 1);
    REQUIRE(formula["O"] == 3);
    REQUIRE(formula["Z"] == -2);

    formula = parseChemicalFormula("Fe+++");
    REQUIRE(formula.size() == 2);
    REQUIRE(formula["Fe"] == 1);
    REQUIRE(formula["Z"] == 3);

    formula = parseChemicalFormula("Fe+3");
    REQUIRE(formula.size() == 2);
    REQUIRE(formula["Fe"] == 1);
    REQUIRE(formula["Z"] == 3);

    formula = parseChemicalFormula("(CaMg)(CO3)2");
    REQUIRE(formula.size() == 4);
    REQUIRE(formula["C"] == 2);
    REQUIRE(formula["Ca"] == 1);
    REQUIRE(formula["Mg"] == 1);
    REQUIRE(formula["O"] == 6);

    formula = parseChemicalFormula("CH3COOH");
    REQUIRE(formula.size() == 3);
    REQUIRE(formula["C"] == 2);
    REQUIRE(formula["H"] == 4);
    REQUIRE(formula["O"] == 2);

    formula = parseChemicalFormula("Al2.5Si0.5O4.75");
    REQUIRE(formula.size() == 3);
    REQUIRE(formula["Al"] == 2.5);
    REQUIRE(formula["Si"] == 0.5);
    REQUIRE(formula["O"] == 4.75);

    formula = parseChemicalFormula("Fe4Al18Si7.5O48H4");
    REQUIRE(formula.size() == 5);
    REQUIRE(formula["Fe"] == 4);
    REQUIRE(formula["Al"] == 18);
    REQUIRE(formula["Si"] == 7.5);
    REQUIRE(formula["O"] == 48);
    REQUIRE(formula["H"] == 4);

    formula = parseChemicalFormula("Mg4Al18Si7.5O48H4");
    REQUIRE(formula.size() == 5);
    REQUIRE(formula["Mg"] == 4);
    REQUIRE(formula["Al"] == 18);
    REQUIRE(formula["Si"] == 7.5);
    REQUIRE(formula["O"] == 48);
    REQUIRE(formula["H"] == 4);

    formula = parseChemicalFormula("Mn4Al18Si7.5O48H4");
    REQUIRE(formula.size() == 5);
    REQUIRE(formula["Mn"] == 4);
    REQUIRE(formula["Al"] == 18);
    REQUIRE(formula["Si"] == 7.5);
    REQUIRE(formula["O"] == 48);
    REQUIRE(formula["H"] == 4);

    formula = parseChemicalFormula("Ca0.5Al1Si2O6");
    REQUIRE(formula.size() == 4);
    REQUIRE(formula["Ca"] == 0.5);
    REQUIRE(formula["Al"] == 1);
    REQUIRE(formula["Si"] == 2);
    REQUIRE(formula["O"] == 6);

    formula = parseChemicalFormula("K0.5Fe5Al2Si8O30.5H12.5");
    REQUIRE(formula.size() == 6);
    REQUIRE(formula["K"] == 0.5);
    REQUIRE(formula["Fe"] == 5);
    REQUIRE(formula["Al"] == 2);
    REQUIRE(formula["Si"] == 8);
    REQUIRE(formula["O"] == 30.5);
    REQUIRE(formula["H"] == 12.5);

    formula = parseChemicalFormula("K0.5Mg5Al2Si8O30.5H12.5");
    REQUIRE(formula.size() == 6);
    REQUIRE(formula["K"] == 0.5);
    REQUIRE(formula["Mg"] == 5);
    REQUIRE(formula["Al"] == 2);
    REQUIRE(formula["Si"] == 8);
    REQUIRE(formula["O"] == 30.5);
    REQUIRE(formula["H"] == 12.5);

    formula = parseChemicalFormula("Mg3.5Al9Si1.5O20");
    REQUIRE(formula.size() == 4);
    REQUIRE(formula["Mg"] == 3.5);
    REQUIRE(formula["Al"] == 9);
    REQUIRE(formula["Si"] == 1.5);
    REQUIRE(formula["O"] == 20);

    formula = parseChemicalFormula("Fe3.5Al9Si1.5O20");
    REQUIRE(formula.size() == 4);
    REQUIRE(formula["Fe"] == 3.5);
    REQUIRE(formula["Al"] == 9);
    REQUIRE(formula["Si"] == 1.5);
    REQUIRE(formula["O"] == 20);

    formula = parseChemicalFormula("Fe0.875S1");
    REQUIRE(formula.size() == 2);
    REQUIRE(formula["Fe"] == 0.875);
    REQUIRE(formula["S"] == 1);
}

TEST_CASE("Testing parseElectricCharge method", "[ChemicalFormula]")
{
    REQUIRE(parseElectricCharge("H+") == +1);
    REQUIRE(parseElectricCharge("H+1") == +1);
    REQUIRE(parseElectricCharge("H(+)") == +1);
    REQUIRE(parseElectricCharge("H+1(aq)") == +1);
    REQUIRE(parseElectricCharge("H+(aq)") == +1);
    REQUIRE(parseElectricCharge("H(+)(aq)") == +1);

    REQUIRE(parseElectricCharge("Ca++") == +2);
    REQUIRE(parseElectricCharge("Ca+2") == +2);
    REQUIRE(parseElectricCharge("Ca(2+)") == +2);
    REQUIRE(parseElectricCharge("Ca++(aq)") == +2);
    REQUIRE(parseElectricCharge("Ca+2(aq)") == +2);
    REQUIRE(parseElectricCharge("Ca(2+)(aq)") == +2);

    REQUIRE(parseElectricCharge("CO3--") == -2);
    REQUIRE(parseElectricCharge("CO3-2") == -2);
    REQUIRE(parseElectricCharge("CO3(2-)") == -2);
    REQUIRE(parseElectricCharge("CO3--(aq)") == -2);
    REQUIRE(parseElectricCharge("CO3-2(aq)") == -2);
    REQUIRE(parseElectricCharge("CO3(2-)(aq)") == -2);
}
