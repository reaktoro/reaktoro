// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/SpeciesList.hpp>
using namespace Reaktoro;

TEST_CASE("Testing SpeciesList", "[SpeciesList]")
{
    SpeciesList species;
    SpeciesList filtered;

    // Test constructor SpeciesList(formulas)
    species = SpeciesList("H2O H+ OH- H2 O2 Na+ Cl- NaCl CO2 HCO3- CO3-2 CH4");

    REQUIRE( species.size() == 12 );

    REQUIRE( species[0].name()  == "H2O"   );
    REQUIRE( species[1].name()  == "H+"    );
    REQUIRE( species[2].name()  == "OH-"   );
    REQUIRE( species[3].name()  == "H2"    );
    REQUIRE( species[4].name()  == "O2"    );
    REQUIRE( species[5].name()  == "Na+"   );
    REQUIRE( species[6].name()  == "Cl-"   );
    REQUIRE( species[7].name()  == "NaCl"  );
    REQUIRE( species[8].name()  == "CO2"   );
    REQUIRE( species[9].name()  == "HCO3-" );
    REQUIRE( species[10].name() == "CO3-2" );
    REQUIRE( species[11].name() == "CH4"   );

    // Test constructor SpeciesList(vector<Species>)
    species = SpeciesList({
        Species("H2O(aq)"   ).withTags({ "aqueous", "neutral", "solvent" }),
        Species("H+(aq)"    ).withTags({ "aqueous", "charged", "cation"})  ,
        Species("OH-(aq)"   ).withTags({ "aqueous", "charged", "anion" })  ,
        Species("H2(aq)"    ).withTags({ "aqueous", "neutral" })           ,
        Species("O2(aq)"    ).withTags({ "aqueous", "neutral" })           ,
        Species("Na+(aq)"   ).withTags({ "aqueous", "charged", "cation"})  ,
        Species("Cl-(aq)"   ).withTags({ "aqueous", "charged", "anion" })  ,
        Species("NaCl(aq)"  ).withTags({ "aqueous", "neutral" })           ,
        Species("CO2(aq)"   ).withTags({ "aqueous", "neutral" })           ,
        Species("HCO3-(aq)" ).withTags({ "aqueous", "charged", "anion" })  ,
        Species("CO3-2(aq)" ).withTags({ "aqueous", "charged", "anion" })  ,
        Species("CH4(aq)"   ).withTags({ "aqueous", "neutral" })           ,
        Species("H2O(g)"    ).withTags({ "gaseous" })                      ,
        Species("CO2(g)"    ).withTags({ "gaseous" })                      ,
        Species("CH4(g)"    ).withTags({ "gaseous" })                      ,
    });

    REQUIRE(species.size() == 15);

    // Test method SpeciesList::indexWithName
    REQUIRE( species.indexWithName("H2O(aq)")   == 0);
    REQUIRE( species.indexWithName("H+(aq)")    == 1);
    REQUIRE( species.indexWithName("OH-(aq)")   == 2);
    REQUIRE( species.indexWithName("H2(aq)")    == 3);
    REQUIRE( species.indexWithName("O2(aq)")    == 4);
    REQUIRE( species.indexWithName("Na+(aq)")   == 5);
    REQUIRE( species.indexWithName("Cl-(aq)")   == 6);
    REQUIRE( species.indexWithName("NaCl(aq)")  == 7);
    REQUIRE( species.indexWithName("CO2(aq)")   == 8);
    REQUIRE( species.indexWithName("HCO3-(aq)") == 9);
    REQUIRE( species.indexWithName("CO3-2(aq)") == 10);
    REQUIRE( species.indexWithName("CH4(aq)")   == 11);
    REQUIRE( species.indexWithName("H2O(g)")    == 12);
    REQUIRE( species.indexWithName("CO2(g)")    == 13);
    REQUIRE( species.indexWithName("CH4(g)")    == 14);

    // Test method SpeciesList::indexWithFormula
    REQUIRE( species.indexWithFormula("H2O")   == 0);
    REQUIRE( species.indexWithFormula("H+")    == 1);
    REQUIRE( species.indexWithFormula("OH-")   == 2);
    REQUIRE( species.indexWithFormula("H2")    == 3);
    REQUIRE( species.indexWithFormula("O2")    == 4);
    REQUIRE( species.indexWithFormula("Na+")   == 5);
    REQUIRE( species.indexWithFormula("Cl-")   == 6);
    REQUIRE( species.indexWithFormula("NaCl")  == 7);
    REQUIRE( species.indexWithFormula("CO2")   == 8);
    REQUIRE( species.indexWithFormula("HCO3-") == 9);
    REQUIRE( species.indexWithFormula("CO3-2") == 10);
    REQUIRE( species.indexWithFormula("CH4")   == 11);

    //-------------------------------------------------------------------------
    // Test method SpeciesList::withNames
    //-------------------------------------------------------------------------
    filtered = species.withNames("H+(aq) OH-(aq) H2O(aq) CO2(g)");

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("H+(aq)" ) == 0 );
    REQUIRE( filtered.indexWithName("OH-(aq)") == 1 );
    REQUIRE( filtered.indexWithName("H2O(aq)") == 2 );
    REQUIRE( filtered.indexWithName("CO2(g)")  == 3 );

    REQUIRE( species.withNames("").size() == 0 );

    REQUIRE_THROWS( species.withNames("ABC") );
    REQUIRE_THROWS( species.withNames("H+(aq) ABC") );

    //-------------------------------------------------------------------------
    // Test method SpeciesList::withFormulas
    //-------------------------------------------------------------------------
    filtered = species.withFormulas("CO2 HCO3- CO3-2 CH4");

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("CO2(aq)"  ) == 0 );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") == 1 );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") == 2 );
    REQUIRE( filtered.indexWithName("CH4(aq)")   == 3 );

    REQUIRE( species.withFormulas("").size() == 0 );

    REQUIRE_THROWS( species.withFormulas("AaBbCc") );
    REQUIRE_THROWS( species.withFormulas("H2O AaBbCc") );

    //-------------------------------------------------------------------------
    // Test method SpeciesList::withSubstances
    //-------------------------------------------------------------------------
    filtered = species.withSubstances("H+ OH- H2O CO2");

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("H+(aq)" ) == 0 );
    REQUIRE( filtered.indexWithName("OH-(aq)") == 1 );
    REQUIRE( filtered.indexWithName("H2O(aq)") == 2 );
    REQUIRE( filtered.indexWithName("CO2(aq)") == 3 );

    REQUIRE( species.withSubstances("").size() == 0 );

    REQUIRE_THROWS( species.withSubstances("ABC") );
    REQUIRE_THROWS( species.withSubstances("H+(aq) ABC") );

    //-------------------------------------------------------------------------
    // Test method SpeciesList::withTag
    //-------------------------------------------------------------------------
    filtered = species.withTag("charged");

    REQUIRE( filtered.size() == 6 );
    REQUIRE( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Na+(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // Test method SpeciesList::withTags
    //-------------------------------------------------------------------------
    filtered = species.withTags({"cation", "charged"});

    REQUIRE( filtered.size() == 2 );
    REQUIRE( filtered.indexWithName("H+(aq)" ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Na+(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // Test method SpeciesList::withTags
    //-------------------------------------------------------------------------
    filtered = species.withTags({"anion", "charged"});

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // Test method SpeciesList::withTag
    //-------------------------------------------------------------------------
    filtered = species.withTag("aqueous");

    REQUIRE( filtered.size() == 12 );
    REQUIRE( filtered.indexWithName("H2O(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("H2(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("O2(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Na+(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("NaCl(aq)" ) < filtered.size() );
    REQUIRE( filtered.indexWithName("CO2(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CH4(aq)"  ) < filtered.size() );

    // Test method SpeciesList::withoutTag
    filtered = species.withoutTag("aqueous");

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered.indexWithName("H2O(g)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO2(g)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CH4(g)") < filtered.size() );

    // Test method SpeciesList::withElements
    filtered = species.withElements("H O");

    REQUIRE( filtered.size() == 6 );
    REQUIRE( filtered.indexWithName("H2O(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H+(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("O2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H2O(g)") < filtered.size() );

    // Test method SpeciesList::withElements
    filtered = species.withElements("H O C");

    REQUIRE( filtered.size() == 12 );
    REQUIRE( filtered.indexWithName("H2O(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("H2(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("O2(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("CO2(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CH4(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("H2O(g)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("CO2(g)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("CH4(g)"   ) < filtered.size() );

    REQUIRE( species.withElements("").size() == 0);
    REQUIRE( species.withElements("Aa Bb Cc Dd Ee").size() == 0 );
    REQUIRE( species.withElements("H O Aa Bb Cc Dd Ee").size() == 6 ); // { H2O(aq), H+(aq), OH-(aq), H2(aq), O2(aq), H2O(g) }

    // Test method SpeciesList::withElementsOf
    filtered = species.withElementsOf({"H2O", "NaCl"});

    REQUIRE( filtered.size() == 9 );
    REQUIRE( filtered.indexWithName("H2O(aq)")   < filtered.size() );
    REQUIRE( filtered.indexWithName("H+(aq)")    < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)")   < filtered.size() );
    REQUIRE( filtered.indexWithName("H2(aq)")    < filtered.size() );
    REQUIRE( filtered.indexWithName("O2(aq)")    < filtered.size() );
    REQUIRE( filtered.indexWithName("Na+(aq)")   < filtered.size() );
    REQUIRE( filtered.indexWithName("Cl-(aq)")   < filtered.size() );
    REQUIRE( filtered.indexWithName("NaCl(aq)")  < filtered.size() );
    REQUIRE( filtered.indexWithName("H2O(g)")    < filtered.size() );

    filtered = filtered.withElementsOf("H2O Na+ Cl-");

    REQUIRE( filtered.size() == 9 );
    REQUIRE( filtered.indexWithName("H2O(aq)")   < filtered.size() );
    REQUIRE( filtered.indexWithName("H+(aq)")    < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)")   < filtered.size() );
    REQUIRE( filtered.indexWithName("H2(aq)")    < filtered.size() );
    REQUIRE( filtered.indexWithName("O2(aq)")    < filtered.size() );
    REQUIRE( filtered.indexWithName("Na+(aq)")   < filtered.size() );
    REQUIRE( filtered.indexWithName("Cl-(aq)")   < filtered.size() );
    REQUIRE( filtered.indexWithName("NaCl(aq)")  < filtered.size() );
    REQUIRE( filtered.indexWithName("H2O(g)")    < filtered.size() );

    filtered = species.withElementsOf("HOCNaCl");

    REQUIRE( filtered.size() == species.size() );
    REQUIRE( filtered.indexWithName("H2O(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("H2(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("O2(aq)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Na+(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("NaCl(aq)" ) < filtered.size() );
    REQUIRE( filtered.indexWithName("CO2(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CH4(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("H2O(g)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("CO2(g)"   ) < filtered.size() );
    REQUIRE( filtered.indexWithName("CH4(g)"   ) < filtered.size() );

    // Test method SpeciesList::append
    species.append(Species("CaCO3(calcite)"));

    REQUIRE( species.indexWithName("CaCO3(calcite)") < species.size() );
    REQUIRE( species.indexWithFormula("CaCO3") < species.size() );
}
