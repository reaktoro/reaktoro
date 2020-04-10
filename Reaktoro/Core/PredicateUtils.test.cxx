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
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/PredicateUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing PredicateUtils", "[PredicateUtils]")
{
    std::vector<Species> species;
    std::vector<Element> elements;

    species = {
        Species("H2O(aq)"  ).withTags({ "aqueous", "neutral", "solvent" }),
        Species("H+(aq)"   ).withTags({ "aqueous", "charged", "cation"})  ,
        Species("OH-(aq)"  ).withTags({ "aqueous", "charged", "anion" })  ,
        Species("H2(aq)"   ).withTags({ "aqueous", "neutral" })           ,
        Species("O2(aq)"   ).withTags({ "aqueous", "neutral" })           ,
        Species("Na+(aq)"  ).withTags({ "aqueous", "charged", "cation"})  ,
        Species("Cl-(aq)"  ).withTags({ "aqueous", "charged", "anion" })  ,
        Species("NaCl(aq)" ).withTags({ "aqueous", "neutral" })           ,
        Species("CO2(aq)"  ).withTags({ "aqueous", "neutral" })           ,
        Species("HCO3-(aq)").withTags({ "aqueous", "charged", "anion" })  ,
        Species("CO3-2(aq)").withTags({ "aqueous", "charged", "anion" })  ,
        Species("CH4(aq)"  ).withTags({ "aqueous", "neutral" })           ,
        Species("H2O(g)"   ).withTags({ "gaseous" })                      ,
        Species("CO2(g)"   ).withTags({ "gaseous" })                      ,
        Species("CH4(g)"   ).withTags({ "gaseous" })                      ,
    };

    elements = {
        Element("H"),
        Element("O"),
        Element("C"),
        Element("Na"),
        Element("Cl"),
        Element("Ca"),
        Element("Mg"),
        Element("Si"),
        Element("Fe"),
    };

    // Test predicate function withSymbol
    REQUIRE( indexfn(elements, withSymbol("H"))  == 0 );
    REQUIRE( indexfn(elements, withSymbol("O"))  == 1 );
    REQUIRE( indexfn(elements, withSymbol("C"))  == 2 );
    REQUIRE( indexfn(elements, withSymbol("Na")) == 3 );
    REQUIRE( indexfn(elements, withSymbol("Cl")) == 4 );
    REQUIRE( indexfn(elements, withSymbol("Ca")) == 5 );
    REQUIRE( indexfn(elements, withSymbol("Mg")) == 6 );
    REQUIRE( indexfn(elements, withSymbol("Si")) == 7 );
    REQUIRE( indexfn(elements, withSymbol("Fe")) == 8 );

    // Test predicate function withName
    REQUIRE( indexfn(species, withName("H2O(aq)"))   == 0);
    REQUIRE( indexfn(species, withName("H+(aq)"))    == 1);
    REQUIRE( indexfn(species, withName("OH-(aq)"))   == 2);
    REQUIRE( indexfn(species, withName("H2(aq)"))    == 3);
    REQUIRE( indexfn(species, withName("O2(aq)"))    == 4);
    REQUIRE( indexfn(species, withName("Na+(aq)"))   == 5);
    REQUIRE( indexfn(species, withName("Cl-(aq)"))   == 6);
    REQUIRE( indexfn(species, withName("NaCl(aq)"))  == 7);
    REQUIRE( indexfn(species, withName("CO2(aq)"))   == 8);
    REQUIRE( indexfn(species, withName("HCO3-(aq)")) == 9);
    REQUIRE( indexfn(species, withName("CO3-2(aq)")) == 10);
    REQUIRE( indexfn(species, withName("CH4(aq)"))   == 11);
    REQUIRE( indexfn(species, withName("H2O(g)"))    == 12);
    REQUIRE( indexfn(species, withName("CO2(g)"))    == 13);
    REQUIRE( indexfn(species, withName("CH4(g)"))    == 14);

    // Test predicate function withFormula
    REQUIRE( indexfn(species, withFormula("H2O"))    == 0);
    REQUIRE( indexfn(species, withFormula("H+"))     == 1);
    REQUIRE( indexfn(species, withFormula("OH-"))    == 2);
    REQUIRE( indexfn(species, withFormula("H2"))     == 3);
    REQUIRE( indexfn(species, withFormula("O2"))     == 4);
    REQUIRE( indexfn(species, withFormula("Na+"))    == 5);
    REQUIRE( indexfn(species, withFormula("Cl-"))    == 6);
    REQUIRE( indexfn(species, withFormula("NaCl"))   == 7);
    REQUIRE( indexfn(species, withFormula("CO2"))    == 8);
    REQUIRE( indexfn(species, withFormula("HCO3-"))  == 9);
    REQUIRE( indexfn(species, withFormula("CO3-2"))  == 10);
    REQUIRE( indexfn(species, withFormula("CH4"))    == 11);

    std::vector<Species> filtered;

    // Test predicate function withNames
    filtered = filter(species, withNames({"H+(aq)", "ABC", "OH-(aq)", "XYZ", "H2O(aq)"}));

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered[0].name() == "H+(aq)"  );
    REQUIRE( filtered[1].name() == "OH-(aq)" );
    REQUIRE( filtered[2].name() == "H2O(aq)" );

    // Test predicate function withFormulas
    filtered = filter(species, withFormulas({"H+", "ABC", "OH-", "XYZ", "H2O"}));

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered[0].name() == "H+(aq)"  );
    REQUIRE( filtered[1].name() == "OH-(aq)" );
    REQUIRE( filtered[2].name() == "H2O(aq)" );

    // Test predicate function withFormulas
    filtered = filter(species, withFormulas({"CO2", "AaBbCc", "HCO3-", "Xx2Yy3", "CO3-2"}));

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered[0].name() == "CO2(aq)"   );
    REQUIRE( filtered[1].name() == "HCO3-(aq)" );
    REQUIRE( filtered[2].name() == "CO3-2(aq)" );

    // Test predicate function withTag
    filtered = filter(species, withTag("charged"));

    REQUIRE( filtered.size() == 6 );
    REQUIRE( filtered[0].name() == "H+(aq)"    );
    REQUIRE( filtered[1].name() == "OH-(aq)"   );
    REQUIRE( filtered[2].name() == "Na+(aq)"   );
    REQUIRE( filtered[3].name() == "Cl-(aq)"   );
    REQUIRE( filtered[4].name() == "HCO3-(aq)" );
    REQUIRE( filtered[5].name() == "CO3-2(aq)" );

    // Test predicate function withTags
    filtered = filter(species, withTags({"cation", "charged"}));

    REQUIRE( filtered.size() == 2 );
    REQUIRE( filtered[0].name() == "H+(aq)"  );
    REQUIRE( filtered[1].name() == "Na+(aq)" );

    // Test predicate function withTags
    filtered = filter(species, withTags({"anion", "charged"}));

    REQUIRE( filtered.size() == 4 );
    REQUIRE( filtered[0].name() == "OH-(aq)"   );
    REQUIRE( filtered[1].name() == "Cl-(aq)"   );
    REQUIRE( filtered[2].name() == "HCO3-(aq)" );
    REQUIRE( filtered[3].name() == "CO3-2(aq)" );

    // Test predicate function withTag
    filtered = filter(species, withTag("aqueous"));

    REQUIRE( filtered.size() == 12 );
    REQUIRE( filtered[0].name()  == "H2O(aq)"   );
    REQUIRE( filtered[1].name()  == "H+(aq)"    );
    REQUIRE( filtered[2].name()  == "OH-(aq)"   );
    REQUIRE( filtered[3].name()  == "H2(aq)"    );
    REQUIRE( filtered[4].name()  == "O2(aq)"    );
    REQUIRE( filtered[5].name()  == "Na+(aq)"   );
    REQUIRE( filtered[6].name()  == "Cl-(aq)"   );
    REQUIRE( filtered[7].name()  == "NaCl(aq)"  );
    REQUIRE( filtered[8].name()  == "CO2(aq)"   );
    REQUIRE( filtered[9].name()  == "HCO3-(aq)" );
    REQUIRE( filtered[10].name() == "CO3-2(aq)" );
    REQUIRE( filtered[11].name() == "CH4(aq)"   );

    // Test predicate function withoutTag
    filtered = filter(species, withoutTag("aqueous"));

    REQUIRE( filtered.size() == 3 );
    REQUIRE( filtered[0].name() == "H2O(g)" );
    REQUIRE( filtered[1].name() == "CO2(g)" );
    REQUIRE( filtered[2].name() == "CH4(g)" );

    // Test predicate function withElements
    filtered = filter(species, withElements({"H", "O"}));

    REQUIRE( filtered.size() == 6 );
    REQUIRE( filtered[0].name() == "H2O(aq)" );
    REQUIRE( filtered[1].name() == "H+(aq)" );
    REQUIRE( filtered[2].name() == "OH-(aq)" );
    REQUIRE( filtered[3].name() == "H2(aq)" );
    REQUIRE( filtered[4].name() == "O2(aq)" );
    REQUIRE( filtered[5].name() == "H2O(g)" );

    // Test predicate function withElements
    filtered = filter(species, withElements({"H", "O", "C"}));

    REQUIRE( filtered.size() == 12 );
    REQUIRE( filtered[0].name()  == "H2O(aq)"   );
    REQUIRE( filtered[1].name()  == "H+(aq)"    );
    REQUIRE( filtered[2].name()  == "OH-(aq)"   );
    REQUIRE( filtered[3].name()  == "H2(aq)"    );
    REQUIRE( filtered[4].name()  == "O2(aq)"    );
    REQUIRE( filtered[5].name()  == "CO2(aq)"   );
    REQUIRE( filtered[6].name()  == "HCO3-(aq)" );
    REQUIRE( filtered[7].name()  == "CO3-2(aq)" );
    REQUIRE( filtered[8].name()  == "CH4(aq)"   );
    REQUIRE( filtered[9].name()  == "H2O(g)"    );
    REQUIRE( filtered[10].name() == "CO2(g)"    );
    REQUIRE( filtered[11].name() == "CH4(g)"    );

    REQUIRE( filter(species, withElements({})).size() == 0 );
    REQUIRE( filter(species, withElements({"Aa", "Bb", "Cc", "Dd", "Ee"})).size() == 0 );
    REQUIRE( filter(species, withElements({"H", "O", "Aa", "Bb", "Cc", "Dd", "Ee"})).size() == 6 ); // { H2O(aq), H+(aq), OH-(aq), H2(aq), O2(aq), H2O(g) }

    // Test predicate function withElementsOf
    filtered = filter(species, withElementsOf({"H2O", "NaCl"}));

    REQUIRE( filtered.size() == 5 );
    REQUIRE( filtered[0].name() == "H2O(aq)"   );
    REQUIRE( filtered[1].name() == "H2(aq)"    );
    REQUIRE( filtered[2].name() == "O2(aq)"    );
    REQUIRE( filtered[3].name() == "NaCl(aq)"  );
    REQUIRE( filtered[4].name() == "H2O(g)"    );

    filtered = filter(species, withElementsOf({"H2O", "Na+", "Cl-"}));

    REQUIRE( filtered.size() == 9 );
    REQUIRE( filtered[0].name() == "H2O(aq)"   );
    REQUIRE( filtered[1].name() == "H+(aq)"    );
    REQUIRE( filtered[2].name() == "OH-(aq)"   );
    REQUIRE( filtered[3].name() == "H2(aq)"    );
    REQUIRE( filtered[4].name() == "O2(aq)"    );
    REQUIRE( filtered[5].name() == "Na+(aq)"   );
    REQUIRE( filtered[6].name() == "Cl-(aq)"   );
    REQUIRE( filtered[7].name() == "NaCl(aq)"  );
    REQUIRE( filtered[8].name() == "H2O(g)"    );

    filtered = filter(species, withElementsOf({"HOCNaClZ"}));

    REQUIRE( filtered.size() == species.size() );
    REQUIRE( filtered[0].name()  == "H2O(aq)"   );
    REQUIRE( filtered[1].name()  == "H+(aq)"    );
    REQUIRE( filtered[2].name()  == "OH-(aq)"   );
    REQUIRE( filtered[3].name()  == "H2(aq)"    );
    REQUIRE( filtered[4].name()  == "O2(aq)"    );
    REQUIRE( filtered[5].name()  == "Na+(aq)"   );
    REQUIRE( filtered[6].name()  == "Cl-(aq)"   );
    REQUIRE( filtered[7].name()  == "NaCl(aq)"  );
    REQUIRE( filtered[8].name()  == "CO2(aq)"   );
    REQUIRE( filtered[9].name()  == "HCO3-(aq)" );
    REQUIRE( filtered[10].name() == "CO3-2(aq)" );
    REQUIRE( filtered[11].name() == "CH4(aq)"   );
    REQUIRE( filtered[12].name() == "H2O(g)"    );
    REQUIRE( filtered[13].name() == "CO2(g)"    );
    REQUIRE( filtered[14].name() == "CH4(g)"    );
}
