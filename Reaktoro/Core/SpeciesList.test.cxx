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
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
using namespace Reaktoro;

TEST_CASE("Testing SpeciesList", "[SpeciesList]")
{
    SpeciesList specieslist;
    SpeciesList filtered;

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: SpeciesList(formulas)
    //-------------------------------------------------------------------------
    specieslist = SpeciesList("H2O H+ OH- H2 O2 Na+ Cl- NaCl CO2 HCO3- CO3-2 CH4");

    CHECK( specieslist.size() == 12 );

    CHECK( specieslist[0].name()  == "H2O"   );
    CHECK( specieslist[1].name()  == "H+"    );
    CHECK( specieslist[2].name()  == "OH-"   );
    CHECK( specieslist[3].name()  == "H2"    );
    CHECK( specieslist[4].name()  == "O2"    );
    CHECK( specieslist[5].name()  == "Na+"   );
    CHECK( specieslist[6].name()  == "Cl-"   );
    CHECK( specieslist[7].name()  == "NaCl"  );
    CHECK( specieslist[8].name()  == "CO2"   );
    CHECK( specieslist[9].name()  == "HCO3-" );
    CHECK( specieslist[10].name() == "CO3-2" );
    CHECK( specieslist[11].name() == "CH4"   );

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: SpeciesList(Vec<Species>)
    //-------------------------------------------------------------------------
    specieslist = SpeciesList({
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

    CHECK(specieslist.size() == 15);

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::elements
    //-------------------------------------------------------------------------
    const auto elements = specieslist.elements();

    CHECK( elements.size() == 5 );
    CHECK( elements[0].symbol() == "H" );
    CHECK( elements[1].symbol() == "C" );
    CHECK( elements[2].symbol() == "O" );
    CHECK( elements[3].symbol() == "Na" );
    CHECK( elements[4].symbol() == "Cl" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::find
    //-------------------------------------------------------------------------
    CHECK( specieslist.find("H2O(aq)")   == 0);
    CHECK( specieslist.find("H+(aq)")    == 1);
    CHECK( specieslist.find("OH-(aq)")   == 2);
    CHECK( specieslist.find("H2(aq)")    == 3);
    CHECK( specieslist.find("O2(aq)")    == 4);
    CHECK( specieslist.find("Na+(aq)")   == 5);
    CHECK( specieslist.find("Cl-(aq)")   == 6);
    CHECK( specieslist.find("NaCl(aq)")  == 7);
    CHECK( specieslist.find("CO2(aq)")   == 8);
    CHECK( specieslist.find("HCO3-(aq)") == 9);
    CHECK( specieslist.find("CO3-2(aq)") == 10);
    CHECK( specieslist.find("CH4(aq)")   == 11);
    CHECK( specieslist.find("H2O(g)")    == 12);
    CHECK( specieslist.find("CO2(g)")    == 13);
    CHECK( specieslist.find("CH4(g)")    == 14);

    CHECK( specieslist.find("XYZ(g)") >= specieslist.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::findWithName
    //-------------------------------------------------------------------------
    CHECK( specieslist.findWithName("H2O(aq)")   == 0);
    CHECK( specieslist.findWithName("H+(aq)")    == 1);
    CHECK( specieslist.findWithName("OH-(aq)")   == 2);
    CHECK( specieslist.findWithName("H2(aq)")    == 3);
    CHECK( specieslist.findWithName("O2(aq)")    == 4);
    CHECK( specieslist.findWithName("Na+(aq)")   == 5);
    CHECK( specieslist.findWithName("Cl-(aq)")   == 6);
    CHECK( specieslist.findWithName("NaCl(aq)")  == 7);
    CHECK( specieslist.findWithName("CO2(aq)")   == 8);
    CHECK( specieslist.findWithName("HCO3-(aq)") == 9);
    CHECK( specieslist.findWithName("CO3-2(aq)") == 10);
    CHECK( specieslist.findWithName("CH4(aq)")   == 11);
    CHECK( specieslist.findWithName("H2O(g)")    == 12);
    CHECK( specieslist.findWithName("CO2(g)")    == 13);
    CHECK( specieslist.findWithName("CH4(g)")    == 14);

    CHECK( specieslist.findWithName("XYZ(g)") >= specieslist.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::findWithFormula
    //-------------------------------------------------------------------------
    CHECK( specieslist.findWithFormula("H2O")   == 0);
    CHECK( specieslist.findWithFormula("H+")    == 1);
    CHECK( specieslist.findWithFormula("OH-")   == 2);
    CHECK( specieslist.findWithFormula("H2")    == 3);
    CHECK( specieslist.findWithFormula("O2")    == 4);
    CHECK( specieslist.findWithFormula("Na+")   == 5);
    CHECK( specieslist.findWithFormula("Cl-")   == 6);
    CHECK( specieslist.findWithFormula("NaCl")  == 7);
    CHECK( specieslist.findWithFormula("CO2")   == 8);
    CHECK( specieslist.findWithFormula("HCO3-") == 9);
    CHECK( specieslist.findWithFormula("CO3-2") == 10);
    CHECK( specieslist.findWithFormula("CH4")   == 11);

    CHECK( specieslist.findWithFormula("XYZ") >= specieslist.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::index
    //-------------------------------------------------------------------------
    CHECK( specieslist.index("H2O(aq)")   == 0);
    CHECK( specieslist.index("H+(aq)")    == 1);
    CHECK( specieslist.index("OH-(aq)")   == 2);
    CHECK( specieslist.index("H2(aq)")    == 3);
    CHECK( specieslist.index("O2(aq)")    == 4);
    CHECK( specieslist.index("Na+(aq)")   == 5);
    CHECK( specieslist.index("Cl-(aq)")   == 6);
    CHECK( specieslist.index("NaCl(aq)")  == 7);
    CHECK( specieslist.index("CO2(aq)")   == 8);
    CHECK( specieslist.index("HCO3-(aq)") == 9);
    CHECK( specieslist.index("CO3-2(aq)") == 10);
    CHECK( specieslist.index("CH4(aq)")   == 11);
    CHECK( specieslist.index("H2O(g)")    == 12);
    CHECK( specieslist.index("CO2(g)")    == 13);
    CHECK( specieslist.index("CH4(g)")    == 14);

    REQUIRE_THROWS( specieslist.index("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::indexWithName
    //-------------------------------------------------------------------------
    CHECK( specieslist.indexWithName("H2O(aq)")   == 0);
    CHECK( specieslist.indexWithName("H+(aq)")    == 1);
    CHECK( specieslist.indexWithName("OH-(aq)")   == 2);
    CHECK( specieslist.indexWithName("H2(aq)")    == 3);
    CHECK( specieslist.indexWithName("O2(aq)")    == 4);
    CHECK( specieslist.indexWithName("Na+(aq)")   == 5);
    CHECK( specieslist.indexWithName("Cl-(aq)")   == 6);
    CHECK( specieslist.indexWithName("NaCl(aq)")  == 7);
    CHECK( specieslist.indexWithName("CO2(aq)")   == 8);
    CHECK( specieslist.indexWithName("HCO3-(aq)") == 9);
    CHECK( specieslist.indexWithName("CO3-2(aq)") == 10);
    CHECK( specieslist.indexWithName("CH4(aq)")   == 11);
    CHECK( specieslist.indexWithName("H2O(g)")    == 12);
    CHECK( specieslist.indexWithName("CO2(g)")    == 13);
    CHECK( specieslist.indexWithName("CH4(g)")    == 14);

    REQUIRE_THROWS( specieslist.indexWithName("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::indexWithFormula
    //-------------------------------------------------------------------------
    CHECK( specieslist.indexWithFormula("H2O")   == 0);
    CHECK( specieslist.indexWithFormula("H+")    == 1);
    CHECK( specieslist.indexWithFormula("OH-")   == 2);
    CHECK( specieslist.indexWithFormula("H2")    == 3);
    CHECK( specieslist.indexWithFormula("O2")    == 4);
    CHECK( specieslist.indexWithFormula("Na+")   == 5);
    CHECK( specieslist.indexWithFormula("Cl-")   == 6);
    CHECK( specieslist.indexWithFormula("NaCl")  == 7);
    CHECK( specieslist.indexWithFormula("CO2")   == 8);
    CHECK( specieslist.indexWithFormula("HCO3-") == 9);
    CHECK( specieslist.indexWithFormula("CO3-2") == 10);
    CHECK( specieslist.indexWithFormula("CH4")   == 11);

    REQUIRE_THROWS( specieslist.indexWithFormula("XYZ") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::get
    //-------------------------------------------------------------------------
    CHECK( specieslist.get("H2O(aq)")  .name() == "H2O(aq)"   );
    CHECK( specieslist.get("H+(aq)")   .name() == "H+(aq)"    );
    CHECK( specieslist.get("OH-(aq)")  .name() == "OH-(aq)"   );
    CHECK( specieslist.get("H2(aq)")   .name() == "H2(aq)"    );
    CHECK( specieslist.get("O2(aq)")   .name() == "O2(aq)"    );
    CHECK( specieslist.get("Na+(aq)")  .name() == "Na+(aq)"   );
    CHECK( specieslist.get("Cl-(aq)")  .name() == "Cl-(aq)"   );
    CHECK( specieslist.get("NaCl(aq)") .name() == "NaCl(aq)"  );
    CHECK( specieslist.get("CO2(aq)")  .name() == "CO2(aq)"   );
    CHECK( specieslist.get("HCO3-(aq)").name() == "HCO3-(aq)" );
    CHECK( specieslist.get("CO3-2(aq)").name() == "CO3-2(aq)" );
    CHECK( specieslist.get("CH4(aq)")  .name() == "CH4(aq)"   );
    CHECK( specieslist.get("H2O(g)")   .name() == "H2O(g)"    );
    CHECK( specieslist.get("CO2(g)")   .name() == "CO2(g)"    );
    CHECK( specieslist.get("CH4(g)")   .name() == "CH4(g)"    );

    REQUIRE_THROWS( specieslist.get("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::getWithName
    //-------------------------------------------------------------------------
    CHECK( specieslist.getWithName("H2O(aq)")  .name() == "H2O(aq)"   );
    CHECK( specieslist.getWithName("H+(aq)")   .name() == "H+(aq)"    );
    CHECK( specieslist.getWithName("OH-(aq)")  .name() == "OH-(aq)"   );
    CHECK( specieslist.getWithName("H2(aq)")   .name() == "H2(aq)"    );
    CHECK( specieslist.getWithName("O2(aq)")   .name() == "O2(aq)"    );
    CHECK( specieslist.getWithName("Na+(aq)")  .name() == "Na+(aq)"   );
    CHECK( specieslist.getWithName("Cl-(aq)")  .name() == "Cl-(aq)"   );
    CHECK( specieslist.getWithName("NaCl(aq)") .name() == "NaCl(aq)"  );
    CHECK( specieslist.getWithName("CO2(aq)")  .name() == "CO2(aq)"   );
    CHECK( specieslist.getWithName("HCO3-(aq)").name() == "HCO3-(aq)" );
    CHECK( specieslist.getWithName("CO3-2(aq)").name() == "CO3-2(aq)" );
    CHECK( specieslist.getWithName("CH4(aq)")  .name() == "CH4(aq)"   );
    CHECK( specieslist.getWithName("H2O(g)")   .name() == "H2O(g)"    );
    CHECK( specieslist.getWithName("CO2(g)")   .name() == "CO2(g)"    );
    CHECK( specieslist.getWithName("CH4(g)")   .name() == "CH4(g)"    );

    REQUIRE_THROWS( specieslist.getWithName("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::getWithFormula
    //-------------------------------------------------------------------------
    CHECK( specieslist.getWithFormula("H2O")  .formula() == "H2O"   );
    CHECK( specieslist.getWithFormula("H+")   .formula() == "H+"    );
    CHECK( specieslist.getWithFormula("OH-")  .formula() == "OH-"   );
    CHECK( specieslist.getWithFormula("H2")   .formula() == "H2"    );
    CHECK( specieslist.getWithFormula("O2")   .formula() == "O2"    );
    CHECK( specieslist.getWithFormula("Na+")  .formula() == "Na+"   );
    CHECK( specieslist.getWithFormula("Cl-")  .formula() == "Cl-"   );
    CHECK( specieslist.getWithFormula("NaCl") .formula() == "NaCl"  );
    CHECK( specieslist.getWithFormula("CO2")  .formula() == "CO2"   );
    CHECK( specieslist.getWithFormula("HCO3-").formula() == "HCO3-" );
    CHECK( specieslist.getWithFormula("CO3-2").formula() == "CO3-2" );
    CHECK( specieslist.getWithFormula("CH4")  .formula() == "CH4"   );

    REQUIRE_THROWS( specieslist.getWithFormula("XYZ") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withNames
    //-------------------------------------------------------------------------
    filtered = specieslist.withNames("H+(aq) OH-(aq) H2O(aq) CO2(g)");

    CHECK( filtered.size() == 4 );
    CHECK( filtered.indexWithName("H+(aq)" ) == 0 );
    CHECK( filtered.indexWithName("OH-(aq)") == 1 );
    CHECK( filtered.indexWithName("H2O(aq)") == 2 );
    CHECK( filtered.indexWithName("CO2(g)")  == 3 );

    CHECK( specieslist.withNames("").size() == 0 );

    REQUIRE_THROWS( specieslist.withNames("ABC") );
    REQUIRE_THROWS( specieslist.withNames("H+(aq) ABC") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withFormulas
    //-------------------------------------------------------------------------
    filtered = specieslist.withFormulas("CO2 HCO3- CO3-2 CH4");

    CHECK( filtered.size() == 4 );
    CHECK( filtered.indexWithName("CO2(aq)"  ) == 0 );
    CHECK( filtered.indexWithName("HCO3-(aq)") == 1 );
    CHECK( filtered.indexWithName("CO3-2(aq)") == 2 );
    CHECK( filtered.indexWithName("CH4(aq)")   == 3 );

    CHECK( specieslist.withFormulas("").size() == 0 );

    REQUIRE_THROWS( specieslist.withFormulas("AaBbCc") );
    REQUIRE_THROWS( specieslist.withFormulas("H2O AaBbCc") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withSubstances
    //-------------------------------------------------------------------------
    filtered = specieslist.withSubstances("H+ OH- H2O CO2");

    CHECK( filtered.size() == 4 );
    CHECK( filtered.indexWithName("H+(aq)" ) == 0 );
    CHECK( filtered.indexWithName("OH-(aq)") == 1 );
    CHECK( filtered.indexWithName("H2O(aq)") == 2 );
    CHECK( filtered.indexWithName("CO2(aq)") == 3 );

    CHECK( specieslist.withSubstances("").size() == 0 );

    REQUIRE_THROWS( specieslist.withSubstances("ABC") );
    REQUIRE_THROWS( specieslist.withSubstances("H+(aq) ABC") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withCharge
    //-------------------------------------------------------------------------
    filtered = specieslist.withCharge(0.0);

    CHECK( filtered.size() == 9 );
    CHECK( filtered[0].name() == "H2O(aq)"  );
    CHECK( filtered[1].name() == "H2(aq)"   );
    CHECK( filtered[2].name() == "O2(aq)"   );
    CHECK( filtered[3].name() == "NaCl(aq)" );
    CHECK( filtered[4].name() == "CO2(aq)"  );
    CHECK( filtered[5].name() == "CH4(aq)"  );
    CHECK( filtered[6].name() == "H2O(g)"   );
    CHECK( filtered[7].name() == "CO2(g)"   );
    CHECK( filtered[8].name() == "CH4(g)"   );

    filtered = specieslist.withCharge(-1.0);

    CHECK( filtered.size() == 3 );
    CHECK( filtered[0].name() == "OH-(aq)"   );
    CHECK( filtered[1].name() == "Cl-(aq)"   );
    CHECK( filtered[2].name() == "HCO3-(aq)" );

    filtered = specieslist.withCharge(-2.0);

    CHECK( filtered.size() == 1 );
    CHECK( filtered[0].name() == "CO3-2(aq)" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTag
    //-------------------------------------------------------------------------
    filtered = specieslist.withTag("charged");

    CHECK( filtered.size() == 6 );
    CHECK( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("Na+(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CO3-2(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTag
    //-------------------------------------------------------------------------
    filtered = specieslist.withTag("charged");

    CHECK( filtered.size() == 6 );
    CHECK( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("Na+(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CO3-2(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTags
    //-------------------------------------------------------------------------
    filtered = specieslist.withTags({"cation", "charged"});

    CHECK( filtered.size() == 2 );
    CHECK( filtered.indexWithName("H+(aq)" ) < filtered.size() );
    CHECK( filtered.indexWithName("Na+(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTags
    //-------------------------------------------------------------------------
    filtered = specieslist.withTags({"anion", "charged"});

    CHECK( filtered.size() == 4 );
    CHECK( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CO3-2(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTag
    //-------------------------------------------------------------------------
    filtered = specieslist.withTag("aqueous");

    CHECK( filtered.size() == 12 );
    CHECK( filtered.indexWithName("H2O(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("H2(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("O2(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("Na+(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("NaCl(aq)" ) < filtered.size() );
    CHECK( filtered.indexWithName("CO2(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CO3-2(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CH4(aq)"  ) < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withoutTag
    //-------------------------------------------------------------------------
    filtered = specieslist.withoutTag("aqueous");

    CHECK( filtered.size() == 3 );
    CHECK( filtered.indexWithName("H2O(g)") < filtered.size() );
    CHECK( filtered.indexWithName("CO2(g)") < filtered.size() );
    CHECK( filtered.indexWithName("CH4(g)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withElements
    //-------------------------------------------------------------------------
    filtered = specieslist.withElements("H O");

    CHECK( filtered.size() == 6 );
    CHECK( filtered.indexWithName("H2O(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("H+(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("OH-(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("H2(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("O2(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("H2O(g)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withElements
    //-------------------------------------------------------------------------
    filtered = specieslist.withElements("H O C");

    CHECK( filtered.size() == 12 );
    CHECK( filtered.indexWithName("H2O(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("H2(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("O2(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("CO2(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CO3-2(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CH4(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("H2O(g)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("CO2(g)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("CH4(g)"   ) < filtered.size() );

    CHECK( specieslist.withElements("").size() == 0);
    CHECK( specieslist.withElements("Aa Bb Cc Dd Ee").size() == 0 );
    CHECK( specieslist.withElements("H O Aa Bb Cc Dd Ee").size() == 6 ); // { H2O(aq), H+(aq), OH-(aq), H2(aq), O2(aq), H2O(g) }

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withElementsOf
    //-------------------------------------------------------------------------
    filtered = specieslist.withElementsOf({"H2O", "NaCl"});

    CHECK( filtered.size() == 9 );
    CHECK( filtered.indexWithName("H2O(aq)")   < filtered.size() );
    CHECK( filtered.indexWithName("H+(aq)")    < filtered.size() );
    CHECK( filtered.indexWithName("OH-(aq)")   < filtered.size() );
    CHECK( filtered.indexWithName("H2(aq)")    < filtered.size() );
    CHECK( filtered.indexWithName("O2(aq)")    < filtered.size() );
    CHECK( filtered.indexWithName("Na+(aq)")   < filtered.size() );
    CHECK( filtered.indexWithName("Cl-(aq)")   < filtered.size() );
    CHECK( filtered.indexWithName("NaCl(aq)")  < filtered.size() );
    CHECK( filtered.indexWithName("H2O(g)")    < filtered.size() );

    filtered = filtered.withElementsOf("H2O Na+ Cl-");

    CHECK( filtered.size() == 9 );
    CHECK( filtered.indexWithName("H2O(aq)")   < filtered.size() );
    CHECK( filtered.indexWithName("H+(aq)")    < filtered.size() );
    CHECK( filtered.indexWithName("OH-(aq)")   < filtered.size() );
    CHECK( filtered.indexWithName("H2(aq)")    < filtered.size() );
    CHECK( filtered.indexWithName("O2(aq)")    < filtered.size() );
    CHECK( filtered.indexWithName("Na+(aq)")   < filtered.size() );
    CHECK( filtered.indexWithName("Cl-(aq)")   < filtered.size() );
    CHECK( filtered.indexWithName("NaCl(aq)")  < filtered.size() );
    CHECK( filtered.indexWithName("H2O(g)")    < filtered.size() );

    filtered = specieslist.withElementsOf("HOCNaCl");

    CHECK( filtered.size() == specieslist.size() );
    CHECK( filtered.indexWithName("H2O(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("H+(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("H2(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("O2(aq)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("Na+(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("NaCl(aq)" ) < filtered.size() );
    CHECK( filtered.indexWithName("CO2(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CO3-2(aq)") < filtered.size() );
    CHECK( filtered.indexWithName("CH4(aq)"  ) < filtered.size() );
    CHECK( filtered.indexWithName("H2O(g)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("CO2(g)"   ) < filtered.size() );
    CHECK( filtered.indexWithName("CH4(g)"   ) < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::append
    //-------------------------------------------------------------------------
    specieslist.append(Species("CaCO3(calcite)"));

    CHECK( specieslist.indexWithName("CaCO3(calcite)") < specieslist.size() );
    CHECK( specieslist.indexWithFormula("CaCO3") < specieslist.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::begin|end
    //-------------------------------------------------------------------------
    for(auto [i, species] : enumerate(specieslist))
        REQUIRE( species.name() == specieslist[i].name() );
}
