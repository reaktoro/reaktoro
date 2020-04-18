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

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: SpeciesList(formulas)
    //-------------------------------------------------------------------------
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

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: SpeciesList(Vec<Species>)
    //-------------------------------------------------------------------------
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

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::find
    //-------------------------------------------------------------------------
    REQUIRE( species.find("H2O(aq)")   == 0);
    REQUIRE( species.find("H+(aq)")    == 1);
    REQUIRE( species.find("OH-(aq)")   == 2);
    REQUIRE( species.find("H2(aq)")    == 3);
    REQUIRE( species.find("O2(aq)")    == 4);
    REQUIRE( species.find("Na+(aq)")   == 5);
    REQUIRE( species.find("Cl-(aq)")   == 6);
    REQUIRE( species.find("NaCl(aq)")  == 7);
    REQUIRE( species.find("CO2(aq)")   == 8);
    REQUIRE( species.find("HCO3-(aq)") == 9);
    REQUIRE( species.find("CO3-2(aq)") == 10);
    REQUIRE( species.find("CH4(aq)")   == 11);
    REQUIRE( species.find("H2O(g)")    == 12);
    REQUIRE( species.find("CO2(g)")    == 13);
    REQUIRE( species.find("CH4(g)")    == 14);

    REQUIRE( species.find("XYZ(g)") >= species.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::findWithName
    //-------------------------------------------------------------------------
    REQUIRE( species.findWithName("H2O(aq)")   == 0);
    REQUIRE( species.findWithName("H+(aq)")    == 1);
    REQUIRE( species.findWithName("OH-(aq)")   == 2);
    REQUIRE( species.findWithName("H2(aq)")    == 3);
    REQUIRE( species.findWithName("O2(aq)")    == 4);
    REQUIRE( species.findWithName("Na+(aq)")   == 5);
    REQUIRE( species.findWithName("Cl-(aq)")   == 6);
    REQUIRE( species.findWithName("NaCl(aq)")  == 7);
    REQUIRE( species.findWithName("CO2(aq)")   == 8);
    REQUIRE( species.findWithName("HCO3-(aq)") == 9);
    REQUIRE( species.findWithName("CO3-2(aq)") == 10);
    REQUIRE( species.findWithName("CH4(aq)")   == 11);
    REQUIRE( species.findWithName("H2O(g)")    == 12);
    REQUIRE( species.findWithName("CO2(g)")    == 13);
    REQUIRE( species.findWithName("CH4(g)")    == 14);

    REQUIRE( species.findWithName("XYZ(g)") >= species.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::findWithFormula
    //-------------------------------------------------------------------------
    REQUIRE( species.findWithFormula("H2O")   == 0);
    REQUIRE( species.findWithFormula("H+")    == 1);
    REQUIRE( species.findWithFormula("OH-")   == 2);
    REQUIRE( species.findWithFormula("H2")    == 3);
    REQUIRE( species.findWithFormula("O2")    == 4);
    REQUIRE( species.findWithFormula("Na+")   == 5);
    REQUIRE( species.findWithFormula("Cl-")   == 6);
    REQUIRE( species.findWithFormula("NaCl")  == 7);
    REQUIRE( species.findWithFormula("CO2")   == 8);
    REQUIRE( species.findWithFormula("HCO3-") == 9);
    REQUIRE( species.findWithFormula("CO3-2") == 10);
    REQUIRE( species.findWithFormula("CH4")   == 11);

    REQUIRE( species.findWithFormula("XYZ") >= species.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::index
    //-------------------------------------------------------------------------
    REQUIRE( species.index("H2O(aq)")   == 0);
    REQUIRE( species.index("H+(aq)")    == 1);
    REQUIRE( species.index("OH-(aq)")   == 2);
    REQUIRE( species.index("H2(aq)")    == 3);
    REQUIRE( species.index("O2(aq)")    == 4);
    REQUIRE( species.index("Na+(aq)")   == 5);
    REQUIRE( species.index("Cl-(aq)")   == 6);
    REQUIRE( species.index("NaCl(aq)")  == 7);
    REQUIRE( species.index("CO2(aq)")   == 8);
    REQUIRE( species.index("HCO3-(aq)") == 9);
    REQUIRE( species.index("CO3-2(aq)") == 10);
    REQUIRE( species.index("CH4(aq)")   == 11);
    REQUIRE( species.index("H2O(g)")    == 12);
    REQUIRE( species.index("CO2(g)")    == 13);
    REQUIRE( species.index("CH4(g)")    == 14);

    REQUIRE_THROWS( species.index("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::indexWithName
    //-------------------------------------------------------------------------
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

    REQUIRE_THROWS( species.indexWithName("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::indexWithFormula
    //-------------------------------------------------------------------------
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

    REQUIRE_THROWS( species.indexWithFormula("XYZ") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withNames
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
    // TESTING METHOD: SpeciesList::withFormulas
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
    // TESTING METHOD: SpeciesList::withSubstances
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
    // TESTING METHOD: SpeciesList::withTag
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
    // TESTING METHOD: SpeciesList::withTags
    //-------------------------------------------------------------------------
    filtered = species.withTags({"cation", "charged"});

    REQUIRE( filtered.size() == 2 );
    REQUIRE( filtered.indexWithName("H+(aq)" ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Na+(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTags
    //-------------------------------------------------------------------------
    filtered = species.withTags({"anion", "charged"});

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTag
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

    // TESTING METHOD: SpeciesList::withoutTag
    filtered = species.withoutTag("aqueous");

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered.indexWithName("H2O(g)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO2(g)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CH4(g)") < filtered.size() );

    // TESTING METHOD: SpeciesList::withElements
    filtered = species.withElements("H O");

    REQUIRE( filtered.size() == 6 );
    REQUIRE( filtered.indexWithName("H2O(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H+(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("O2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H2O(g)") < filtered.size() );

    // TESTING METHOD: SpeciesList::withElements
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

    // TESTING METHOD: SpeciesList::withElementsOf
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

    // TESTING METHOD: SpeciesList::append
    species.append(Species("CaCO3(calcite)"));

    REQUIRE( species.indexWithName("CaCO3(calcite)") < species.size() );
    REQUIRE( species.indexWithFormula("CaCO3") < species.size() );
}
