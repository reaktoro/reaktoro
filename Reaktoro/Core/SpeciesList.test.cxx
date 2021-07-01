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

    REQUIRE( specieslist.size() == 12 );

    REQUIRE( specieslist[0].name()  == "H2O"   );
    REQUIRE( specieslist[1].name()  == "H+"    );
    REQUIRE( specieslist[2].name()  == "OH-"   );
    REQUIRE( specieslist[3].name()  == "H2"    );
    REQUIRE( specieslist[4].name()  == "O2"    );
    REQUIRE( specieslist[5].name()  == "Na+"   );
    REQUIRE( specieslist[6].name()  == "Cl-"   );
    REQUIRE( specieslist[7].name()  == "NaCl"  );
    REQUIRE( specieslist[8].name()  == "CO2"   );
    REQUIRE( specieslist[9].name()  == "HCO3-" );
    REQUIRE( specieslist[10].name() == "CO3-2" );
    REQUIRE( specieslist[11].name() == "CH4"   );

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

    REQUIRE(specieslist.size() == 15);

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::find
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.find("H2O(aq)")   == 0);
    REQUIRE( specieslist.find("H+(aq)")    == 1);
    REQUIRE( specieslist.find("OH-(aq)")   == 2);
    REQUIRE( specieslist.find("H2(aq)")    == 3);
    REQUIRE( specieslist.find("O2(aq)")    == 4);
    REQUIRE( specieslist.find("Na+(aq)")   == 5);
    REQUIRE( specieslist.find("Cl-(aq)")   == 6);
    REQUIRE( specieslist.find("NaCl(aq)")  == 7);
    REQUIRE( specieslist.find("CO2(aq)")   == 8);
    REQUIRE( specieslist.find("HCO3-(aq)") == 9);
    REQUIRE( specieslist.find("CO3-2(aq)") == 10);
    REQUIRE( specieslist.find("CH4(aq)")   == 11);
    REQUIRE( specieslist.find("H2O(g)")    == 12);
    REQUIRE( specieslist.find("CO2(g)")    == 13);
    REQUIRE( specieslist.find("CH4(g)")    == 14);

    REQUIRE( specieslist.find("XYZ(g)") >= specieslist.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::findWithName
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.findWithName("H2O(aq)")   == 0);
    REQUIRE( specieslist.findWithName("H+(aq)")    == 1);
    REQUIRE( specieslist.findWithName("OH-(aq)")   == 2);
    REQUIRE( specieslist.findWithName("H2(aq)")    == 3);
    REQUIRE( specieslist.findWithName("O2(aq)")    == 4);
    REQUIRE( specieslist.findWithName("Na+(aq)")   == 5);
    REQUIRE( specieslist.findWithName("Cl-(aq)")   == 6);
    REQUIRE( specieslist.findWithName("NaCl(aq)")  == 7);
    REQUIRE( specieslist.findWithName("CO2(aq)")   == 8);
    REQUIRE( specieslist.findWithName("HCO3-(aq)") == 9);
    REQUIRE( specieslist.findWithName("CO3-2(aq)") == 10);
    REQUIRE( specieslist.findWithName("CH4(aq)")   == 11);
    REQUIRE( specieslist.findWithName("H2O(g)")    == 12);
    REQUIRE( specieslist.findWithName("CO2(g)")    == 13);
    REQUIRE( specieslist.findWithName("CH4(g)")    == 14);

    REQUIRE( specieslist.findWithName("XYZ(g)") >= specieslist.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::findWithFormula
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.findWithFormula("H2O")   == 0);
    REQUIRE( specieslist.findWithFormula("H+")    == 1);
    REQUIRE( specieslist.findWithFormula("OH-")   == 2);
    REQUIRE( specieslist.findWithFormula("H2")    == 3);
    REQUIRE( specieslist.findWithFormula("O2")    == 4);
    REQUIRE( specieslist.findWithFormula("Na+")   == 5);
    REQUIRE( specieslist.findWithFormula("Cl-")   == 6);
    REQUIRE( specieslist.findWithFormula("NaCl")  == 7);
    REQUIRE( specieslist.findWithFormula("CO2")   == 8);
    REQUIRE( specieslist.findWithFormula("HCO3-") == 9);
    REQUIRE( specieslist.findWithFormula("CO3-2") == 10);
    REQUIRE( specieslist.findWithFormula("CH4")   == 11);

    REQUIRE( specieslist.findWithFormula("XYZ") >= specieslist.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::index
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.index("H2O(aq)")   == 0);
    REQUIRE( specieslist.index("H+(aq)")    == 1);
    REQUIRE( specieslist.index("OH-(aq)")   == 2);
    REQUIRE( specieslist.index("H2(aq)")    == 3);
    REQUIRE( specieslist.index("O2(aq)")    == 4);
    REQUIRE( specieslist.index("Na+(aq)")   == 5);
    REQUIRE( specieslist.index("Cl-(aq)")   == 6);
    REQUIRE( specieslist.index("NaCl(aq)")  == 7);
    REQUIRE( specieslist.index("CO2(aq)")   == 8);
    REQUIRE( specieslist.index("HCO3-(aq)") == 9);
    REQUIRE( specieslist.index("CO3-2(aq)") == 10);
    REQUIRE( specieslist.index("CH4(aq)")   == 11);
    REQUIRE( specieslist.index("H2O(g)")    == 12);
    REQUIRE( specieslist.index("CO2(g)")    == 13);
    REQUIRE( specieslist.index("CH4(g)")    == 14);

    REQUIRE_THROWS( specieslist.index("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::indexWithName
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.indexWithName("H2O(aq)")   == 0);
    REQUIRE( specieslist.indexWithName("H+(aq)")    == 1);
    REQUIRE( specieslist.indexWithName("OH-(aq)")   == 2);
    REQUIRE( specieslist.indexWithName("H2(aq)")    == 3);
    REQUIRE( specieslist.indexWithName("O2(aq)")    == 4);
    REQUIRE( specieslist.indexWithName("Na+(aq)")   == 5);
    REQUIRE( specieslist.indexWithName("Cl-(aq)")   == 6);
    REQUIRE( specieslist.indexWithName("NaCl(aq)")  == 7);
    REQUIRE( specieslist.indexWithName("CO2(aq)")   == 8);
    REQUIRE( specieslist.indexWithName("HCO3-(aq)") == 9);
    REQUIRE( specieslist.indexWithName("CO3-2(aq)") == 10);
    REQUIRE( specieslist.indexWithName("CH4(aq)")   == 11);
    REQUIRE( specieslist.indexWithName("H2O(g)")    == 12);
    REQUIRE( specieslist.indexWithName("CO2(g)")    == 13);
    REQUIRE( specieslist.indexWithName("CH4(g)")    == 14);

    REQUIRE_THROWS( specieslist.indexWithName("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::indexWithFormula
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.indexWithFormula("H2O")   == 0);
    REQUIRE( specieslist.indexWithFormula("H+")    == 1);
    REQUIRE( specieslist.indexWithFormula("OH-")   == 2);
    REQUIRE( specieslist.indexWithFormula("H2")    == 3);
    REQUIRE( specieslist.indexWithFormula("O2")    == 4);
    REQUIRE( specieslist.indexWithFormula("Na+")   == 5);
    REQUIRE( specieslist.indexWithFormula("Cl-")   == 6);
    REQUIRE( specieslist.indexWithFormula("NaCl")  == 7);
    REQUIRE( specieslist.indexWithFormula("CO2")   == 8);
    REQUIRE( specieslist.indexWithFormula("HCO3-") == 9);
    REQUIRE( specieslist.indexWithFormula("CO3-2") == 10);
    REQUIRE( specieslist.indexWithFormula("CH4")   == 11);

    REQUIRE_THROWS( specieslist.indexWithFormula("XYZ") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::get
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.get("H2O(aq)")  .name() == "H2O(aq)"   );
    REQUIRE( specieslist.get("H+(aq)")   .name() == "H+(aq)"    );
    REQUIRE( specieslist.get("OH-(aq)")  .name() == "OH-(aq)"   );
    REQUIRE( specieslist.get("H2(aq)")   .name() == "H2(aq)"    );
    REQUIRE( specieslist.get("O2(aq)")   .name() == "O2(aq)"    );
    REQUIRE( specieslist.get("Na+(aq)")  .name() == "Na+(aq)"   );
    REQUIRE( specieslist.get("Cl-(aq)")  .name() == "Cl-(aq)"   );
    REQUIRE( specieslist.get("NaCl(aq)") .name() == "NaCl(aq)"  );
    REQUIRE( specieslist.get("CO2(aq)")  .name() == "CO2(aq)"   );
    REQUIRE( specieslist.get("HCO3-(aq)").name() == "HCO3-(aq)" );
    REQUIRE( specieslist.get("CO3-2(aq)").name() == "CO3-2(aq)" );
    REQUIRE( specieslist.get("CH4(aq)")  .name() == "CH4(aq)"   );
    REQUIRE( specieslist.get("H2O(g)")   .name() == "H2O(g)"    );
    REQUIRE( specieslist.get("CO2(g)")   .name() == "CO2(g)"    );
    REQUIRE( specieslist.get("CH4(g)")   .name() == "CH4(g)"    );

    REQUIRE_THROWS( specieslist.get("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::getWithName
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.getWithName("H2O(aq)")  .name() == "H2O(aq)"   );
    REQUIRE( specieslist.getWithName("H+(aq)")   .name() == "H+(aq)"    );
    REQUIRE( specieslist.getWithName("OH-(aq)")  .name() == "OH-(aq)"   );
    REQUIRE( specieslist.getWithName("H2(aq)")   .name() == "H2(aq)"    );
    REQUIRE( specieslist.getWithName("O2(aq)")   .name() == "O2(aq)"    );
    REQUIRE( specieslist.getWithName("Na+(aq)")  .name() == "Na+(aq)"   );
    REQUIRE( specieslist.getWithName("Cl-(aq)")  .name() == "Cl-(aq)"   );
    REQUIRE( specieslist.getWithName("NaCl(aq)") .name() == "NaCl(aq)"  );
    REQUIRE( specieslist.getWithName("CO2(aq)")  .name() == "CO2(aq)"   );
    REQUIRE( specieslist.getWithName("HCO3-(aq)").name() == "HCO3-(aq)" );
    REQUIRE( specieslist.getWithName("CO3-2(aq)").name() == "CO3-2(aq)" );
    REQUIRE( specieslist.getWithName("CH4(aq)")  .name() == "CH4(aq)"   );
    REQUIRE( specieslist.getWithName("H2O(g)")   .name() == "H2O(g)"    );
    REQUIRE( specieslist.getWithName("CO2(g)")   .name() == "CO2(g)"    );
    REQUIRE( specieslist.getWithName("CH4(g)")   .name() == "CH4(g)"    );

    REQUIRE_THROWS( specieslist.getWithName("XYZ(g)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::getWithFormula
    //-------------------------------------------------------------------------
    REQUIRE( specieslist.getWithFormula("H2O")  .formula() == "H2O"   );
    REQUIRE( specieslist.getWithFormula("H+")   .formula() == "H+"    );
    REQUIRE( specieslist.getWithFormula("OH-")  .formula() == "OH-"   );
    REQUIRE( specieslist.getWithFormula("H2")   .formula() == "H2"    );
    REQUIRE( specieslist.getWithFormula("O2")   .formula() == "O2"    );
    REQUIRE( specieslist.getWithFormula("Na+")  .formula() == "Na+"   );
    REQUIRE( specieslist.getWithFormula("Cl-")  .formula() == "Cl-"   );
    REQUIRE( specieslist.getWithFormula("NaCl") .formula() == "NaCl"  );
    REQUIRE( specieslist.getWithFormula("CO2")  .formula() == "CO2"   );
    REQUIRE( specieslist.getWithFormula("HCO3-").formula() == "HCO3-" );
    REQUIRE( specieslist.getWithFormula("CO3-2").formula() == "CO3-2" );
    REQUIRE( specieslist.getWithFormula("CH4")  .formula() == "CH4"   );

    REQUIRE_THROWS( specieslist.getWithFormula("XYZ") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withNames
    //-------------------------------------------------------------------------
    filtered = specieslist.withNames("H+(aq) OH-(aq) H2O(aq) CO2(g)");

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("H+(aq)" ) == 0 );
    REQUIRE( filtered.indexWithName("OH-(aq)") == 1 );
    REQUIRE( filtered.indexWithName("H2O(aq)") == 2 );
    REQUIRE( filtered.indexWithName("CO2(g)")  == 3 );

    REQUIRE( specieslist.withNames("").size() == 0 );

    REQUIRE_THROWS( specieslist.withNames("ABC") );
    REQUIRE_THROWS( specieslist.withNames("H+(aq) ABC") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withFormulas
    //-------------------------------------------------------------------------
    filtered = specieslist.withFormulas("CO2 HCO3- CO3-2 CH4");

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("CO2(aq)"  ) == 0 );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") == 1 );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") == 2 );
    REQUIRE( filtered.indexWithName("CH4(aq)")   == 3 );

    REQUIRE( specieslist.withFormulas("").size() == 0 );

    REQUIRE_THROWS( specieslist.withFormulas("AaBbCc") );
    REQUIRE_THROWS( specieslist.withFormulas("H2O AaBbCc") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withSubstances
    //-------------------------------------------------------------------------
    filtered = specieslist.withSubstances("H+ OH- H2O CO2");

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("H+(aq)" ) == 0 );
    REQUIRE( filtered.indexWithName("OH-(aq)") == 1 );
    REQUIRE( filtered.indexWithName("H2O(aq)") == 2 );
    REQUIRE( filtered.indexWithName("CO2(aq)") == 3 );

    REQUIRE( specieslist.withSubstances("").size() == 0 );

    REQUIRE_THROWS( specieslist.withSubstances("ABC") );
    REQUIRE_THROWS( specieslist.withSubstances("H+(aq) ABC") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTag
    //-------------------------------------------------------------------------
    filtered = specieslist.withTag("charged");

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
    filtered = specieslist.withTags({"cation", "charged"});

    REQUIRE( filtered.size() == 2 );
    REQUIRE( filtered.indexWithName("H+(aq)" ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Na+(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTags
    //-------------------------------------------------------------------------
    filtered = specieslist.withTags({"anion", "charged"});

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered.indexWithName("OH-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("Cl-(aq)"  ) < filtered.size() );
    REQUIRE( filtered.indexWithName("HCO3-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO3-2(aq)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withTag
    //-------------------------------------------------------------------------
    filtered = specieslist.withTag("aqueous");

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

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withoutTag
    //-------------------------------------------------------------------------
    filtered = specieslist.withoutTag("aqueous");

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered.indexWithName("H2O(g)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CO2(g)") < filtered.size() );
    REQUIRE( filtered.indexWithName("CH4(g)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withElements
    //-------------------------------------------------------------------------
    filtered = specieslist.withElements("H O");

    REQUIRE( filtered.size() == 6 );
    REQUIRE( filtered.indexWithName("H2O(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H+(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("OH-(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("O2(aq)") < filtered.size() );
    REQUIRE( filtered.indexWithName("H2O(g)") < filtered.size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withElements
    //-------------------------------------------------------------------------
    filtered = specieslist.withElements("H O C");

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

    REQUIRE( specieslist.withElements("").size() == 0);
    REQUIRE( specieslist.withElements("Aa Bb Cc Dd Ee").size() == 0 );
    REQUIRE( specieslist.withElements("H O Aa Bb Cc Dd Ee").size() == 6 ); // { H2O(aq), H+(aq), OH-(aq), H2(aq), O2(aq), H2O(g) }

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::withElementsOf
    //-------------------------------------------------------------------------
    filtered = specieslist.withElementsOf({"H2O", "NaCl"});

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

    filtered = specieslist.withElementsOf("HOCNaCl");

    REQUIRE( filtered.size() == specieslist.size() );
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

    //-------------------------------------------------------------------------
    // TESTING METHOD: SpeciesList::append
    //-------------------------------------------------------------------------
    specieslist.append(Species("CaCO3(calcite)"));

    REQUIRE( specieslist.indexWithName("CaCO3(calcite)") < specieslist.size() );
    REQUIRE( specieslist.indexWithFormula("CaCO3") < specieslist.size() );
}
