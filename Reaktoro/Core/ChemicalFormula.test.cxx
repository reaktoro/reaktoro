// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
    CHECK( formula.str()             == "H2O" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 2 );
    CHECK( formula.coefficient("H")  == 2 );
    CHECK( formula.coefficient("O")  == 1 );
    CHECK( formula.equivalent("H2O(aq)") );
    CHECK( formula.equivalent("H2O(l)") );
    CHECK( formula.equivalent("HHO") );
    CHECK( formula.equivalent("HOH") );
    CHECK( formula.equivalent("OH2") );

    formula = ChemicalFormula("CaCO3");
    CHECK( formula.str()             == "CaCO3" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 3 );
    CHECK( formula.coefficient("C")  == 1 );
    CHECK( formula.coefficient("Ca") == 1 );
    CHECK( formula.coefficient("O")  == 3 );
    CHECK( formula.equivalent("CaCOOO") );
    CHECK( formula.equivalent("CaOOOC") );
    CHECK( formula.equivalent("Ca(CO3)") );
    CHECK( formula.equivalent("Ca(CO3)(aq)") );
    CHECK( formula.equivalent("Ca(CO3)(s)") );

    formula = ChemicalFormula("HCO3-");
    CHECK( formula.str()             == "HCO3-" );
    CHECK( formula.charge()          == -1 );
    CHECK( formula.elements().size() == 3 );
    CHECK( formula.coefficient("C")  == 1 );
    CHECK( formula.coefficient("H")  == 1 );
    CHECK( formula.coefficient("O")  == 3 );
    CHECK( formula.equivalent("HCO3-(aq)") );
    CHECK( formula.equivalent("HCO3[-](aq)") );
    CHECK( formula.equivalent("HCO3[-]") );
    CHECK( formula.equivalent("HCOOO-") );

    formula = ChemicalFormula("H+");
    CHECK( formula.str()             == "H+" );
    CHECK( formula.charge()          == 1 );
    CHECK( formula.elements().size() == 1 );
    CHECK( formula.coefficient("H")  == 1 );
    CHECK( formula.equivalent("H+(aq)") );
    CHECK( formula.equivalent("H[+]") );

    formula = ChemicalFormula("e-");
    CHECK( formula.str()             == "e-" );
    CHECK( formula.charge()          == -1 );
    CHECK( formula.elements().size() == 0 );
    CHECK( formula.equivalent("e-(aq)") );
    CHECK( formula.equivalent("e[-]") );

    formula = ChemicalFormula("Na+");
    CHECK( formula.str()             == "Na+" );
    CHECK( formula.charge()          == 1 );
    CHECK( formula.elements().size() == 1 );
    CHECK( formula.coefficient("Na") == 1 );
    CHECK( formula.equivalent("Na+(aq)") );
    CHECK( formula.equivalent("Na+(pl)") );
    CHECK( formula.equivalent("Na[+](aq)") );
    CHECK( formula.equivalent("Na[+]") );

    formula = ChemicalFormula("Cl-");
    CHECK( formula.str()             == "Cl-" );
    CHECK( formula.charge()          == -1 );
    CHECK( formula.elements().size() == 1 );
    CHECK( formula.coefficient("Cl") == 1 );
    CHECK( formula.equivalent("Cl[-](aq)") );
    CHECK( formula.equivalent("Cl-(aq)") );
    CHECK( formula.equivalent("Cl-(pl)") );

    formula = ChemicalFormula("CO3--");
    CHECK( formula.str()             == "CO3--" );
    CHECK( formula.charge()          == -2 );
    CHECK( formula.elements().size() == 2 );
    CHECK( formula.coefficient("C")  == 1 );
    CHECK( formula.coefficient("O")  == 3 );
    CHECK( formula.equivalent("CO3-2(aq)") );
    CHECK( formula.equivalent("CO3-2") );
    CHECK( formula.equivalent("CO3[2-](aq)") );

    formula = ChemicalFormula("Fe+++");
    CHECK( formula.str()             == "Fe+++" );
    CHECK( formula.charge()          == 3 );
    CHECK( formula.elements().size() == 1 );
    CHECK( formula.coefficient("Fe") == 1 );
    CHECK( formula.equivalent("Fe+++(aq)") );
    CHECK( formula.equivalent("Fe+3(aq)") );
    CHECK( formula.equivalent("Fe[3+](aq)") );
    CHECK( formula.equivalent("Fe+3") );
    CHECK( formula.equivalent("Fe[3+]") );

    formula = ChemicalFormula("(CaMg)(CO3)2");
    CHECK( formula.str()             == "(CaMg)(CO3)2" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 4 );
    CHECK( formula.coefficient("C")  == 2 );
    CHECK( formula.coefficient("Ca") == 1 );
    CHECK( formula.coefficient("Mg") == 1 );
    CHECK( formula.coefficient("O")  == 6 );
    CHECK( formula.equivalent("CaMg(CO3)2(s)") );
    CHECK( formula.equivalent("CaMg(CO3)2") );
    CHECK( formula.equivalent("CaMgCO3CO3(s)") );
    CHECK( formula.equivalent("CaMgCO3CO3") );

    formula = ChemicalFormula("CH3COOH");
    CHECK( formula.str()             == "CH3COOH" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 3 );
    CHECK( formula.coefficient("C")  == 2 );
    CHECK( formula.coefficient("H")  == 4 );
    CHECK( formula.coefficient("O")  == 2 );
    CHECK( formula.equivalent("CH3COOH(aq)") );
    CHECK( formula.equivalent("CHHHCOOH(aq)") );
    CHECK( formula.equivalent("CHHHCOOH") );

    formula = ChemicalFormula("Al2.5Si0.5O4.75");
    CHECK( formula.str()             == "Al2.5Si0.5O4.75" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 3 );
    CHECK( formula.coefficient("Al") == 2.5 );
    CHECK( formula.coefficient("Si") == 0.5 );
    CHECK( formula.coefficient("O")  == 4.75 );
    CHECK( formula.equivalent("Al2.5Si0.5O4.75(s)") );
    CHECK( formula.equivalent("Si0.5Al2.5O4.75(s)") );

    formula = ChemicalFormula("Fe4Al18Si7.5O48H4");
    CHECK( formula.str()             == "Fe4Al18Si7.5O48H4" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 5 );
    CHECK( formula.coefficient("Fe") == 4 );
    CHECK( formula.coefficient("Al") == 18 );
    CHECK( formula.coefficient("Si") == 7.5 );
    CHECK( formula.coefficient("O")  == 48 );
    CHECK( formula.coefficient("H")  == 4 );
    CHECK( formula.equivalent("Fe4Al18Si7.5O48H4(s)") );
    CHECK( formula.equivalent("Al18Fe4H4Si7.5O48(s)") );

    formula = ChemicalFormula("Mg4Al18Si7.5O48H4");
    CHECK( formula.str()             == "Mg4Al18Si7.5O48H4" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 5 );
    CHECK( formula.coefficient("Mg") == 4 );
    CHECK( formula.coefficient("Al") == 18 );
    CHECK( formula.coefficient("Si") == 7.5 );
    CHECK( formula.coefficient("O")  == 48 );
    CHECK( formula.coefficient("H")  == 4 );
    CHECK( formula.equivalent("Mg4Al18Si7.5O48H4(s)") );

    formula = ChemicalFormula("Mn4Al18Si7.5O48H4");
    CHECK( formula.str()             == "Mn4Al18Si7.5O48H4" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 5 );
    CHECK( formula.coefficient("Mn") == 4 );
    CHECK( formula.coefficient("Al") == 18 );
    CHECK( formula.coefficient("Si") == 7.5 );
    CHECK( formula.coefficient("O")  == 48 );
    CHECK( formula.coefficient("H")  == 4 );

    formula = ChemicalFormula("Ca0.5Al1Si2O6");
    CHECK( formula.str()             == "Ca0.5Al1Si2O6" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 4 );
    CHECK( formula.coefficient("Ca") == 0.5 );
    CHECK( formula.coefficient("Al") == 1 );
    CHECK( formula.coefficient("Si") == 2 );
    CHECK( formula.coefficient("O")  == 6 );

    formula = ChemicalFormula("K0.5Fe5Al2Si8O30.5H12.5");
    CHECK( formula.str()             == "K0.5Fe5Al2Si8O30.5H12.5" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 6 );
    CHECK( formula.coefficient("K")  == 0.5 );
    CHECK( formula.coefficient("Fe") == 5 );
    CHECK( formula.coefficient("Al") == 2 );
    CHECK( formula.coefficient("Si") == 8 );
    CHECK( formula.coefficient("O")  == 30.5 );
    CHECK( formula.coefficient("H")  == 12.5 );

    formula = ChemicalFormula("K0.5Mg5Al2Si8O30.5H12.5");
    CHECK( formula.str()             == "K0.5Mg5Al2Si8O30.5H12.5" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 6 );
    CHECK( formula.coefficient("K")  == 0.5 );
    CHECK( formula.coefficient("Mg") == 5 );
    CHECK( formula.coefficient("Al") == 2 );
    CHECK( formula.coefficient("Si") == 8 );
    CHECK( formula.coefficient("O")  == 30.5 );
    CHECK( formula.coefficient("H")  == 12.5 );

    formula = ChemicalFormula("Mg3.5Al9Si1.5O20");
    CHECK( formula.str()             == "Mg3.5Al9Si1.5O20" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 4 );
    CHECK( formula.coefficient("Mg") == 3.5 );
    CHECK( formula.coefficient("Al") == 9 );
    CHECK( formula.coefficient("Si") == 1.5 );
    CHECK( formula.coefficient("O")  == 20 );

    formula = ChemicalFormula("Fe3.5Al9Si1.5O20");
    CHECK( formula.str()             == "Fe3.5Al9Si1.5O20" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 4 );
    CHECK( formula.coefficient("Fe") == 3.5 );
    CHECK( formula.coefficient("Al") == 9 );
    CHECK( formula.coefficient("Si") == 1.5 );
    CHECK( formula.coefficient("O")  == 20 );

    formula = ChemicalFormula("Fe0.875S1");
    CHECK( formula.str()             == "Fe0.875S1" );
    CHECK( formula.charge()          == 0 );
    CHECK( formula.elements().size() == 2 );
    CHECK( formula.coefficient("Fe") == 0.875 );
    CHECK( formula.coefficient("S")  == 1 );

    CHECK(ChemicalFormula::equivalent("Ca++", "Ca+2"));
    CHECK(ChemicalFormula::equivalent("Ca++", "Ca[2+]"));
    CHECK(ChemicalFormula::equivalent("Ca++", "Ca[2+](aq)"));

    CHECK(ChemicalFormula::equivalent("CO3--", "CO3-2"));
    CHECK(ChemicalFormula::equivalent("CO3--", "CO3[2-]"));
    CHECK(ChemicalFormula::equivalent("CO3--", "CO3[2-](aq)"));

    CHECK(ChemicalFormula::equivalent("Fe+++", "Fe+3"));
    CHECK(ChemicalFormula::equivalent("Fe+++", "Fe[3+]"));
    CHECK(ChemicalFormula::equivalent("Fe+++", "Fe[3+](aq)"));

    CHECK(ChemicalFormula::equivalent("H+", "H+1"));
    CHECK(ChemicalFormula::equivalent("H+", "H[+]"));
    CHECK(ChemicalFormula::equivalent("H+", "H[+](aq)"));

    CHECK(ChemicalFormula::equivalent("OH-", "OH-1"));
    CHECK(ChemicalFormula::equivalent("OH-", "OH[-]"));
    CHECK(ChemicalFormula::equivalent("OH-", "OH[-](aq)"));

    CHECK(ChemicalFormula::equivalent("CO2", "CO2(g)"));
    CHECK(ChemicalFormula::equivalent("CO2", "COO"));
}

TEST_CASE("Testing parseChemicalFormula method", "[ChemicalFormula]")
{
    Pairs<String, double> formula;

    formula = parseChemicalFormula("H2O");
    CHECK( formula == Pairs<String, double>{{"H", 2}, {"O", 1}} );

    formula = parseChemicalFormula("CaCO3");
    CHECK( formula == Pairs<String, double>{{"Ca", 1}, {"C", 1}, {"O", 3}} );

    formula = parseChemicalFormula("HCO3-");
    CHECK( formula == Pairs<String, double>{{"H", 1}, {"C", 1}, {"O", 3}} );

    formula = parseChemicalFormula("H+");
    CHECK( formula == Pairs<String, double>{{"H", 1}} );

    formula = parseChemicalFormula("Na+");
    CHECK( formula == Pairs<String, double>{{"Na", 1}} );

    formula = parseChemicalFormula("Cl-");
    CHECK( formula == Pairs<String, double>{{"Cl", 1}} );

    formula = parseChemicalFormula("CO3--");
    CHECK( formula == Pairs<String, double>{{"C", 1}, {"O", 3}} );

    formula = parseChemicalFormula("CO3-2");
    CHECK( formula == Pairs<String, double>{{"C", 1}, {"O", 3}} );

    formula = parseChemicalFormula("Fe+++");
    CHECK( formula == Pairs<String, double>{{"Fe", 1}} );

    formula = parseChemicalFormula("Fe+3");
    CHECK( formula == Pairs<String, double>{{"Fe", 1}} );

    formula = parseChemicalFormula("(CaMg)(CO3)2");
    CHECK( formula == Pairs<String, double>{{"Ca", 1}, {"Mg", 1}, {"C", 2}, {"O", 6}} );

    formula = parseChemicalFormula("CH3COOH");
    CHECK( formula == Pairs<String, double>{{"C", 2}, {"H", 4}, {"O", 2}} );

    formula = parseChemicalFormula("Al2.5Si0.5O4.75");
    CHECK( formula == Pairs<String, double>{{"Al", 2.5}, {"Si", 0.5}, {"O", 4.75}} );

    formula = parseChemicalFormula("Fe4Al18Si7.5O48H4");
    CHECK( formula == Pairs<String, double>{{"Fe", 4}, {"Al", 18}, {"Si", 7.5}, {"O", 48}, {"H", 4}} );

    formula = parseChemicalFormula("Mg4Al18Si7.5O48H4");
    CHECK( formula == Pairs<String, double>{{"Mg", 4}, {"Al", 18}, {"Si", 7.5}, {"O", 48}, {"H", 4}} );

    formula = parseChemicalFormula("Mn4Al18Si7.5O48H4");
    CHECK( formula == Pairs<String, double>{{"Mn", 4}, {"Al", 18}, {"Si", 7.5}, {"O", 48}, {"H", 4}} );

    formula = parseChemicalFormula("Ca0.5Al1Si2O6");
    CHECK( formula == Pairs<String, double>{{"Ca", 0.5}, {"Al", 1}, {"Si", 2}, {"O", 6}} );

    formula = parseChemicalFormula("K0.5Fe5Al2Si8O30.5H12.5");
    CHECK( formula == Pairs<String, double>{{"K", 0.5}, {"Fe", 5}, {"Al", 2}, {"Si", 8}, {"O", 30.5}, {"H", 12.5}} );

    formula = parseChemicalFormula("K0.5Mg5Al2Si8O30.5H12.5");
    CHECK( formula == Pairs<String, double>{{"K", 0.5}, {"Mg", 5}, {"Al", 2}, {"Si", 8}, {"O", 30.5}, {"H", 12.5}} );

    formula = parseChemicalFormula("Mg3.5Al9Si1.5O20");
    CHECK( formula == Pairs<String, double>{{"Mg", 3.5}, {"Al", 9}, {"Si", 1.5}, {"O", 20}} );

    formula = parseChemicalFormula("Fe3.5Al9Si1.5O20");
    CHECK( formula == Pairs<String, double>{{"Fe", 3.5}, {"Al", 9}, {"Si", 1.5}, {"O", 20}} );

    formula = parseChemicalFormula("Fe0.875S1");
    CHECK( formula == Pairs<String, double>{{"Fe", 0.875}, {"S", 1}} );
}

TEST_CASE("Testing parseElectricCharge method", "[ChemicalFormula]")
{
    CHECK( parseElectricCharge("H+") == +1 );
    CHECK( parseElectricCharge("H+1") == +1 );
    CHECK( parseElectricCharge("H[+]") == +1 );
    CHECK( parseElectricCharge("H+1(aq)") == +1 );
    CHECK( parseElectricCharge("H+(aq)") == +1 );
    CHECK( parseElectricCharge("H[+](aq)") == +1 );

    CHECK( parseElectricCharge("Ca++") == +2 );
    CHECK( parseElectricCharge("Ca+2") == +2 );
    CHECK( parseElectricCharge("Ca[2+]") == +2 );
    CHECK( parseElectricCharge("Ca++(aq)") == +2 );
    CHECK( parseElectricCharge("Ca+2(aq)") == +2 );
    CHECK( parseElectricCharge("Ca[2+](aq)") == +2 );

    CHECK( parseElectricCharge("CO3--") == -2 );
    CHECK( parseElectricCharge("CO3-2") == -2 );
    CHECK( parseElectricCharge("CO3[2-]") == -2 );
    CHECK( parseElectricCharge("CO3--(aq)") == -2 );
    CHECK( parseElectricCharge("CO3-2(aq)") == -2 );
    CHECK( parseElectricCharge("CO3[2-](aq)") == -2 );
}
