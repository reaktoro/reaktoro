// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/ParseUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ParseUtils module", "[ParseUtils]")
{
    //-----------------------------------------------------------------
    // Testing method: parseReaction
    //-----------------------------------------------------------------

    // TODO: Rename parseReaction to parseReactionEquation and implement tests for this method.

    //-----------------------------------------------------------------
    // Testing method: parseNumberStringPairs
    //-----------------------------------------------------------------
    Pairs<String, double> pairs;
    pairs = parseNumberStringPairs("2:H 1:O");
    CHECK( pairs == Pairs<String, double>{{"H", 2}, {"O", 1}} );

    pairs = parseNumberStringPairs("1:Ca 2:Cl");
    CHECK( pairs == Pairs<String, double>{{"Ca", 1}, {"Cl", 2}} );

    pairs = parseNumberStringPairs("1:Mg 1:C 3:O");
    CHECK( pairs == Pairs<String, double>{{"Mg", 1}, {"C", 1}, {"O", 3}} );

    pairs = parseNumberStringPairs("1:Ca 1:Mg 2:C 6:O");
    CHECK( pairs == Pairs<String, double>{{"Ca", 1}, {"Mg", 1}, {"C", 2}, {"O", 6}} );

    pairs = parseNumberStringPairs("3:Fe 2:Al 3:Si 12:O");
    CHECK( pairs == Pairs<String, double>{{"Fe", 3}, {"Al", 2}, {"Si", 3}, {"O", 12}} );

    pairs = parseNumberStringPairs("1:Fe 1:Al 1:Fe 3:Si 1:Fe 1:Al 4:O 4:O 4:O");
    CHECK( pairs == Pairs<String, double>{{"Fe", 3}, {"Al", 2}, {"Si", 3}, {"O", 12}} );

    pairs = parseNumberStringPairs("1:Na -1:E");
    CHECK( pairs == Pairs<String, double>{{"Na", 1}, {"E", -1}} );

    pairs = parseNumberStringPairs("1:Na -1:E -3:E 3:E");
    CHECK( pairs == Pairs<String, double>{{"Na", 1}, {"E", -1}} );

    pairs = parseNumberStringPairs("1:Ca -2:E");
    CHECK( pairs == Pairs<String, double>{{"Ca", 1}, {"E", -2}} );

    pairs = parseNumberStringPairs("1:Fe -3:E");
    CHECK( pairs == Pairs<String, double>{{"Fe", 1}, {"E", -3}} );

    pairs = parseNumberStringPairs("1:C 3:O 2:E");
    CHECK( pairs == Pairs<String, double>{{"C", 1}, {"O", 3}, {"E", 2}} );

    pairs = parseNumberStringPairs("1:Ca++ 1:CO3--");
    CHECK( pairs == Pairs<String, double>{{"Ca++", 1}, {"CO3--", 1}} );

    pairs = parseNumberStringPairs("1:Ca++ 1:Mg++ 2:CO3--");
    CHECK( pairs == Pairs<String, double>{{"Ca++", 1}, {"Mg++", 1}, {"CO3--", 2}} );
}

TEST_CASE("Testing parseChemicalFormula method", "[ParseUtils]")
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

    formula = parseChemicalFormula("CaCl2*6H2O");
    CHECK( formula == Pairs<String, double>{{"Ca", 1}, {"Cl", 2}, {"H", 12}, {"O", 6}} );

    formula = parseChemicalFormula("SrCl2*2H2O");
    CHECK( formula == Pairs<String, double>{{"Sr", 1}, {"Cl", 2}, {"H", 4}, {"O", 2}} );

    formula = parseChemicalFormula("(NH4)2SO4*3NH4NO3");
    CHECK( formula == Pairs<String, double>{{"N", 8}, {"H", 20}, {"S", 1}, {"O", 13}} );

    formula = parseChemicalFormula("Na2SO4*(NH4)2SO4*4H2O");
    CHECK( formula == Pairs<String, double>{{"Na", 2}, {"S", 2}, {"O", 12}, {"N", 2}, {"H", 16}} );

    formula = parseChemicalFormula("CaCl2:6H2O");
    CHECK( formula == Pairs<String, double>{{"Ca", 1}, {"Cl", 2}, {"H", 12}, {"O", 6}} );

    formula = parseChemicalFormula("SrCl2:2H2O");
    CHECK( formula == Pairs<String, double>{{"Sr", 1}, {"Cl", 2}, {"H", 4}, {"O", 2}} );

    formula = parseChemicalFormula("(NH4)2SO4:3NH4NO3");
    CHECK( formula == Pairs<String, double>{{"N", 8}, {"H", 20}, {"S", 1}, {"O", 13}} );

    formula = parseChemicalFormula("Na2SO4:(NH4)2SO4*4H2O");
    CHECK( formula == Pairs<String, double>{{"Na", 2}, {"S", 2}, {"O", 12}, {"N", 2}, {"H", 16}} );

    formula = parseChemicalFormula("Na2SO4:(NH4)2SO4:4H2O");
    CHECK( formula == Pairs<String, double>{{"Na", 2}, {"S", 2}, {"O", 12}, {"N", 2}, {"H", 16}} );
}

TEST_CASE("Testing parseElectricCharge method", "[ParseUtils]")
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

TEST_CASE("Testing parseReactionEquation method", "[ParseUtils]")
{
    Pairs<String, double> equation;

    equation = parseReactionEquation("H2O = H+ + OH-");

    REQUIRE( equation.size() == 3            );
    REQUIRE( equation.at(0).first  == "H2O"  );
    REQUIRE( equation.at(1).first  == "H+"   );
    REQUIRE( equation.at(2).first  == "OH-"  );
    REQUIRE( equation.at(0).second == -1.0   );
    REQUIRE( equation.at(1).second ==  1.0   );
    REQUIRE( equation.at(2).second ==  1.0   );

    equation = parseReactionEquation("CaCl2 = Ca++ + 2*Cl-");

    REQUIRE( equation.size() == 3             );
    REQUIRE( equation.at(0).first  == "CaCl2" );
    REQUIRE( equation.at(1).first  == "Ca++"  );
    REQUIRE( equation.at(2).first  == "Cl-"   );
    REQUIRE( equation.at(0).second == -1.0    );
    REQUIRE( equation.at(1).second ==  1.0    );
    REQUIRE( equation.at(2).second ==  2.0    );

    equation = parseReactionEquation("Ca++ + Mg++ + 2*HCO3- = CaMg(CO3)2 + 2*H+");

    REQUIRE( equation.size() == 5                  );
    REQUIRE( equation.at(0).first  == "Ca++"       );
    REQUIRE( equation.at(1).first  == "Mg++"       );
    REQUIRE( equation.at(2).first  == "HCO3-"      );
    REQUIRE( equation.at(3).first  == "CaMg(CO3)2" );
    REQUIRE( equation.at(4).first  == "H+"         );
    REQUIRE( equation.at(0).second == -1.0         );
    REQUIRE( equation.at(1).second == -1.0         );
    REQUIRE( equation.at(2).second == -2.0         );
    REQUIRE( equation.at(3).second ==  1.0         );
    REQUIRE( equation.at(4).second ==  2.0         );

    equation = parseReactionEquation("CaCO3(s)");

    REQUIRE( equation.size() == 1                );
    REQUIRE( equation.at(0).first  == "CaCO3(s)" );
    REQUIRE( equation.at(0).second == -1.0       );

    equation = parseReactionEquation("2*CaCO3(s) + 3*H+");

    REQUIRE( equation.size() == 2                );
    REQUIRE( equation.at(0).first  == "CaCO3(s)" );
    REQUIRE( equation.at(1).first  == "H+"       );
    REQUIRE( equation.at(0).second == -2.0       );
    REQUIRE( equation.at(1).second == -3.0       );
}
