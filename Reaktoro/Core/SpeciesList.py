# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.


from reaktoro import *
import pytest


def testSpeciesList():

    #-------------------------------------------------------------------------
    # TESTING CONSTRUCTOR: SpeciesList(formulas)
    #-------------------------------------------------------------------------
    specieslist = SpeciesList("H2O H+ OH- H2 O2 Na+ Cl- NaCl CO2 HCO3- CO3-2 CH4")

    assert specieslist.size() == 12

    assert specieslist[0].name()  == "H2O"
    assert specieslist[1].name()  == "H+"
    assert specieslist[2].name()  == "OH-"
    assert specieslist[3].name()  == "H2"
    assert specieslist[4].name()  == "O2"
    assert specieslist[5].name()  == "Na+"
    assert specieslist[6].name()  == "Cl-"
    assert specieslist[7].name()  == "NaCl"
    assert specieslist[8].name()  == "CO2"
    assert specieslist[9].name()  == "HCO3-"
    assert specieslist[10].name() == "CO3-2"
    assert specieslist[11].name() == "CH4"

    #-------------------------------------------------------------------------
    # TESTING CONSTRUCTOR: SpeciesList(Vec<Species>)
    #-------------------------------------------------------------------------
    specieslist = SpeciesList([
        Species("H2O(aq)"   ).withTags([ "aqueous", "neutral", "solvent" ]),
        Species("H+(aq)"    ).withTags([ "aqueous", "charged", "cation"])  ,
        Species("OH-(aq)"   ).withTags([ "aqueous", "charged", "anion" ])  ,
        Species("H2(aq)"    ).withTags([ "aqueous", "neutral" ])           ,
        Species("O2(aq)"    ).withTags([ "aqueous", "neutral" ])           ,
        Species("Na+(aq)"   ).withTags([ "aqueous", "charged", "cation"])  ,
        Species("Cl-(aq)"   ).withTags([ "aqueous", "charged", "anion" ])  ,
        Species("NaCl(aq)"  ).withTags([ "aqueous", "neutral" ])           ,
        Species("CO2(aq)"   ).withTags([ "aqueous", "neutral" ])           ,
        Species("HCO3-(aq)" ).withTags([ "aqueous", "charged", "anion" ])  ,
        Species("CO3-2(aq)" ).withTags([ "aqueous", "charged", "anion" ])  ,
        Species("CH4(aq)"   ).withTags([ "aqueous", "neutral" ])           ,
        Species("H2O(g)"    ).withTags([ "gaseous" ])                      ,
        Species("CO2(g)"    ).withTags([ "gaseous" ])                      ,
        Species("CH4(g)"    ).withTags([ "gaseous" ])                      ,
    ])

    assert specieslist.size() == 15

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::find
    #-------------------------------------------------------------------------
    assert specieslist.find("H2O(aq)")   == 0
    assert specieslist.find("H+(aq)")    == 1
    assert specieslist.find("OH-(aq)")   == 2
    assert specieslist.find("H2(aq)")    == 3
    assert specieslist.find("O2(aq)")    == 4
    assert specieslist.find("Na+(aq)")   == 5
    assert specieslist.find("Cl-(aq)")   == 6
    assert specieslist.find("NaCl(aq)")  == 7
    assert specieslist.find("CO2(aq)")   == 8
    assert specieslist.find("HCO3-(aq)") == 9
    assert specieslist.find("CO3-2(aq)") == 10
    assert specieslist.find("CH4(aq)")   == 11
    assert specieslist.find("H2O(g)")    == 12
    assert specieslist.find("CO2(g)")    == 13
    assert specieslist.find("CH4(g)")    == 14

    assert specieslist.find("XYZ(g)") >= specieslist.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::findWithName
    #-------------------------------------------------------------------------
    assert specieslist.findWithName("H2O(aq)")   == 0
    assert specieslist.findWithName("H+(aq)")    == 1
    assert specieslist.findWithName("OH-(aq)")   == 2
    assert specieslist.findWithName("H2(aq)")    == 3
    assert specieslist.findWithName("O2(aq)")    == 4
    assert specieslist.findWithName("Na+(aq)")   == 5
    assert specieslist.findWithName("Cl-(aq)")   == 6
    assert specieslist.findWithName("NaCl(aq)")  == 7
    assert specieslist.findWithName("CO2(aq)")   == 8
    assert specieslist.findWithName("HCO3-(aq)") == 9
    assert specieslist.findWithName("CO3-2(aq)") == 10
    assert specieslist.findWithName("CH4(aq)")   == 11
    assert specieslist.findWithName("H2O(g)")    == 12
    assert specieslist.findWithName("CO2(g)")    == 13
    assert specieslist.findWithName("CH4(g)")    == 14

    assert specieslist.findWithName("XYZ(g)") >= specieslist.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::findWithFormula
    #-------------------------------------------------------------------------
    assert specieslist.findWithFormula("H2O")   == 0
    assert specieslist.findWithFormula("H+")    == 1
    assert specieslist.findWithFormula("OH-")   == 2
    assert specieslist.findWithFormula("H2")    == 3
    assert specieslist.findWithFormula("O2")    == 4
    assert specieslist.findWithFormula("Na+")   == 5
    assert specieslist.findWithFormula("Cl-")   == 6
    assert specieslist.findWithFormula("NaCl")  == 7
    assert specieslist.findWithFormula("CO2")   == 8
    assert specieslist.findWithFormula("HCO3-") == 9
    assert specieslist.findWithFormula("CO3-2") == 10
    assert specieslist.findWithFormula("CH4")   == 11

    assert specieslist.findWithFormula("XYZ") >= specieslist.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::index
    #-------------------------------------------------------------------------
    assert specieslist.index("H2O(aq)")   == 0
    assert specieslist.index("H+(aq)")    == 1
    assert specieslist.index("OH-(aq)")   == 2
    assert specieslist.index("H2(aq)")    == 3
    assert specieslist.index("O2(aq)")    == 4
    assert specieslist.index("Na+(aq)")   == 5
    assert specieslist.index("Cl-(aq)")   == 6
    assert specieslist.index("NaCl(aq)")  == 7
    assert specieslist.index("CO2(aq)")   == 8
    assert specieslist.index("HCO3-(aq)") == 9
    assert specieslist.index("CO3-2(aq)") == 10
    assert specieslist.index("CH4(aq)")   == 11
    assert specieslist.index("H2O(g)")    == 12
    assert specieslist.index("CO2(g)")    == 13
    assert specieslist.index("CH4(g)")    == 14

    with pytest.raises(Exception):
        specieslist.index("XYZ(g)")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::indexWithName
    #-------------------------------------------------------------------------
    assert specieslist.indexWithName("H2O(aq)")   == 0
    assert specieslist.indexWithName("H+(aq)")    == 1
    assert specieslist.indexWithName("OH-(aq)")   == 2
    assert specieslist.indexWithName("H2(aq)")    == 3
    assert specieslist.indexWithName("O2(aq)")    == 4
    assert specieslist.indexWithName("Na+(aq)")   == 5
    assert specieslist.indexWithName("Cl-(aq)")   == 6
    assert specieslist.indexWithName("NaCl(aq)")  == 7
    assert specieslist.indexWithName("CO2(aq)")   == 8
    assert specieslist.indexWithName("HCO3-(aq)") == 9
    assert specieslist.indexWithName("CO3-2(aq)") == 10
    assert specieslist.indexWithName("CH4(aq)")   == 11
    assert specieslist.indexWithName("H2O(g)")    == 12
    assert specieslist.indexWithName("CO2(g)")    == 13
    assert specieslist.indexWithName("CH4(g)")    == 14

    with pytest.raises(Exception):
        specieslist.indexWithName("XYZ(g)")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::indexWithFormula
    #-------------------------------------------------------------------------
    assert specieslist.indexWithFormula("H2O")   == 0
    assert specieslist.indexWithFormula("H+")    == 1
    assert specieslist.indexWithFormula("OH-")   == 2
    assert specieslist.indexWithFormula("H2")    == 3
    assert specieslist.indexWithFormula("O2")    == 4
    assert specieslist.indexWithFormula("Na+")   == 5
    assert specieslist.indexWithFormula("Cl-")   == 6
    assert specieslist.indexWithFormula("NaCl")  == 7
    assert specieslist.indexWithFormula("CO2")   == 8
    assert specieslist.indexWithFormula("HCO3-") == 9
    assert specieslist.indexWithFormula("CO3-2") == 10
    assert specieslist.indexWithFormula("CH4")   == 11

    with pytest.raises(Exception):
        specieslist.indexWithFormula("XYZ")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::get
    #-------------------------------------------------------------------------
    assert specieslist.get("H2O(aq)")  .name() == "H2O(aq)"
    assert specieslist.get("H+(aq)")   .name() == "H+(aq)"
    assert specieslist.get("OH-(aq)")  .name() == "OH-(aq)"
    assert specieslist.get("H2(aq)")   .name() == "H2(aq)"
    assert specieslist.get("O2(aq)")   .name() == "O2(aq)"
    assert specieslist.get("Na+(aq)")  .name() == "Na+(aq)"
    assert specieslist.get("Cl-(aq)")  .name() == "Cl-(aq)"
    assert specieslist.get("NaCl(aq)") .name() == "NaCl(aq)"
    assert specieslist.get("CO2(aq)")  .name() == "CO2(aq)"
    assert specieslist.get("HCO3-(aq)").name() == "HCO3-(aq)"
    assert specieslist.get("CO3-2(aq)").name() == "CO3-2(aq)"
    assert specieslist.get("CH4(aq)")  .name() == "CH4(aq)"
    assert specieslist.get("H2O(g)")   .name() == "H2O(g)"
    assert specieslist.get("CO2(g)")   .name() == "CO2(g)"
    assert specieslist.get("CH4(g)")   .name() == "CH4(g)"

    with pytest.raises(Exception):
        specieslist.get("XYZ(g)")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::getWithName
    #-------------------------------------------------------------------------
    assert specieslist.getWithName("H2O(aq)")  .name() == "H2O(aq)"
    assert specieslist.getWithName("H+(aq)")   .name() == "H+(aq)"
    assert specieslist.getWithName("OH-(aq)")  .name() == "OH-(aq)"
    assert specieslist.getWithName("H2(aq)")   .name() == "H2(aq)"
    assert specieslist.getWithName("O2(aq)")   .name() == "O2(aq)"
    assert specieslist.getWithName("Na+(aq)")  .name() == "Na+(aq)"
    assert specieslist.getWithName("Cl-(aq)")  .name() == "Cl-(aq)"
    assert specieslist.getWithName("NaCl(aq)") .name() == "NaCl(aq)"
    assert specieslist.getWithName("CO2(aq)")  .name() == "CO2(aq)"
    assert specieslist.getWithName("HCO3-(aq)").name() == "HCO3-(aq)"
    assert specieslist.getWithName("CO3-2(aq)").name() == "CO3-2(aq)"
    assert specieslist.getWithName("CH4(aq)")  .name() == "CH4(aq)"
    assert specieslist.getWithName("H2O(g)")   .name() == "H2O(g)"
    assert specieslist.getWithName("CO2(g)")   .name() == "CO2(g)"
    assert specieslist.getWithName("CH4(g)")   .name() == "CH4(g)"

    with pytest.raises(Exception):
        specieslist.getWithName("XYZ(g)")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::getWithFormula
    #-------------------------------------------------------------------------
    assert specieslist.getWithFormula("H2O")  .formula() == "H2O"
    assert specieslist.getWithFormula("H+")   .formula() == "H+"
    assert specieslist.getWithFormula("OH-")  .formula() == "OH-"
    assert specieslist.getWithFormula("H2")   .formula() == "H2"
    assert specieslist.getWithFormula("O2")   .formula() == "O2"
    assert specieslist.getWithFormula("Na+")  .formula() == "Na+"
    assert specieslist.getWithFormula("Cl-")  .formula() == "Cl-"
    assert specieslist.getWithFormula("NaCl") .formula() == "NaCl"
    assert specieslist.getWithFormula("CO2")  .formula() == "CO2"
    assert specieslist.getWithFormula("HCO3-").formula() == "HCO3-"
    assert specieslist.getWithFormula("CO3-2").formula() == "CO3-2"
    assert specieslist.getWithFormula("CH4")  .formula() == "CH4"

    with pytest.raises(Exception):
        specieslist.getWithFormula("XYZ")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withNames
    #-------------------------------------------------------------------------
    filtered = specieslist.withNames("H+(aq) OH-(aq) H2O(aq) CO2(g)")

    assert filtered.size() == 4
    assert filtered.indexWithName("H+(aq)" ) == 0
    assert filtered.indexWithName("OH-(aq)") == 1
    assert filtered.indexWithName("H2O(aq)") == 2
    assert filtered.indexWithName("CO2(g)")  == 3

    assert specieslist.withNames("").size() == 0

    with pytest.raises(Exception):
        specieslist.withNames("ABC")
    with pytest.raises(Exception):
        specieslist.withNames("H+(aq) ABC")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withFormulas
    #-------------------------------------------------------------------------
    filtered = specieslist.withFormulas("CO2 HCO3- CO3-2 CH4")

    assert filtered.size() == 4
    assert filtered.indexWithName("CO2(aq)"  ) == 0
    assert filtered.indexWithName("HCO3-(aq)") == 1
    assert filtered.indexWithName("CO3-2(aq)") == 2
    assert filtered.indexWithName("CH4(aq)")   == 3

    assert specieslist.withFormulas("").size() == 0

    with pytest.raises(Exception):
        specieslist.withFormulas("AaBbCc")
    with pytest.raises(Exception):
        specieslist.withFormulas("H2O AaBbCc")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withSubstances
    #-------------------------------------------------------------------------
    filtered = specieslist.withSubstances("H+ OH- H2O CO2")

    assert filtered.size() == 4
    assert filtered.indexWithName("H+(aq)" ) == 0
    assert filtered.indexWithName("OH-(aq)") == 1
    assert filtered.indexWithName("H2O(aq)") == 2
    assert filtered.indexWithName("CO2(aq)") == 3

    assert specieslist.withSubstances("").size() == 0

    with pytest.raises(Exception):
        specieslist.withSubstances("ABC")
    with pytest.raises(Exception):
        specieslist.withSubstances("H+(aq) ABC")

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withTag
    #-------------------------------------------------------------------------
    filtered = specieslist.withTag("charged")

    assert filtered.size() == 6
    assert filtered.indexWithName("H+(aq)"   ) < filtered.size()
    assert filtered.indexWithName("OH-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("Na+(aq)"  ) < filtered.size()
    assert filtered.indexWithName("Cl-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("HCO3-(aq)") < filtered.size()
    assert filtered.indexWithName("CO3-2(aq)") < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withTags
    #-------------------------------------------------------------------------
    filtered = specieslist.withTags(["cation", "charged"])

    assert filtered.size() == 2
    assert filtered.indexWithName("H+(aq)" ) < filtered.size()
    assert filtered.indexWithName("Na+(aq)") < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withTags
    #-------------------------------------------------------------------------
    filtered = specieslist.withTags(["anion", "charged"])

    assert filtered.size() == 4
    assert filtered.indexWithName("OH-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("Cl-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("HCO3-(aq)") < filtered.size()
    assert filtered.indexWithName("CO3-2(aq)") < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withTag
    #-------------------------------------------------------------------------
    filtered = specieslist.withTag("aqueous")

    assert filtered.size() == 12
    assert filtered.indexWithName("H2O(aq)"  ) < filtered.size()
    assert filtered.indexWithName("H+(aq)"   ) < filtered.size()
    assert filtered.indexWithName("OH-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("H2(aq)"   ) < filtered.size()
    assert filtered.indexWithName("O2(aq)"   ) < filtered.size()
    assert filtered.indexWithName("Na+(aq)"  ) < filtered.size()
    assert filtered.indexWithName("Cl-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("NaCl(aq)" ) < filtered.size()
    assert filtered.indexWithName("CO2(aq)"  ) < filtered.size()
    assert filtered.indexWithName("HCO3-(aq)") < filtered.size()
    assert filtered.indexWithName("CO3-2(aq)") < filtered.size()
    assert filtered.indexWithName("CH4(aq)"  ) < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withoutTag
    #-------------------------------------------------------------------------
    filtered = specieslist.withoutTag("aqueous")

    assert filtered.size() == 3
    assert filtered.indexWithName("H2O(g)") < filtered.size()
    assert filtered.indexWithName("CO2(g)") < filtered.size()
    assert filtered.indexWithName("CH4(g)") < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withElements
    #-------------------------------------------------------------------------
    filtered = specieslist.withElements("H O")

    assert filtered.size() == 6
    assert filtered.indexWithName("H2O(aq)") < filtered.size()
    assert filtered.indexWithName("H+(aq)") < filtered.size()
    assert filtered.indexWithName("OH-(aq)") < filtered.size()
    assert filtered.indexWithName("H2(aq)") < filtered.size()
    assert filtered.indexWithName("O2(aq)") < filtered.size()
    assert filtered.indexWithName("H2O(g)") < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withElements
    #-------------------------------------------------------------------------
    filtered = specieslist.withElements("H O C")

    assert filtered.size() == 12
    assert filtered.indexWithName("H2O(aq)"  ) < filtered.size()
    assert filtered.indexWithName("H+(aq)"   ) < filtered.size()
    assert filtered.indexWithName("OH-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("H2(aq)"   ) < filtered.size()
    assert filtered.indexWithName("O2(aq)"   ) < filtered.size()
    assert filtered.indexWithName("CO2(aq)"  ) < filtered.size()
    assert filtered.indexWithName("HCO3-(aq)") < filtered.size()
    assert filtered.indexWithName("CO3-2(aq)") < filtered.size()
    assert filtered.indexWithName("CH4(aq)"  ) < filtered.size()
    assert filtered.indexWithName("H2O(g)"   ) < filtered.size()
    assert filtered.indexWithName("CO2(g)"   ) < filtered.size()
    assert filtered.indexWithName("CH4(g)"   ) < filtered.size()

    assert specieslist.withElements("").size() == 0
    assert specieslist.withElements("Aa Bb Cc Dd Ee").size() == 0
    assert specieslist.withElements("H O Aa Bb Cc Dd Ee").size() == 6  # { H2O(aq), H+(aq), OH-(aq), H2(aq), O2(aq), H2O(g) }

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::withElementsOf
    #-------------------------------------------------------------------------
    filtered = specieslist.withElementsOf(["H2O", "NaCl"])

    assert filtered.size() == 9
    assert filtered.indexWithName("H2O(aq)")   < filtered.size()
    assert filtered.indexWithName("H+(aq)")    < filtered.size()
    assert filtered.indexWithName("OH-(aq)")   < filtered.size()
    assert filtered.indexWithName("H2(aq)")    < filtered.size()
    assert filtered.indexWithName("O2(aq)")    < filtered.size()
    assert filtered.indexWithName("Na+(aq)")   < filtered.size()
    assert filtered.indexWithName("Cl-(aq)")   < filtered.size()
    assert filtered.indexWithName("NaCl(aq)")  < filtered.size()
    assert filtered.indexWithName("H2O(g)")    < filtered.size()

    filtered = filtered.withElementsOf("H2O Na+ Cl-")

    assert filtered.size() == 9
    assert filtered.indexWithName("H2O(aq)")   < filtered.size()
    assert filtered.indexWithName("H+(aq)")    < filtered.size()
    assert filtered.indexWithName("OH-(aq)")   < filtered.size()
    assert filtered.indexWithName("H2(aq)")    < filtered.size()
    assert filtered.indexWithName("O2(aq)")    < filtered.size()
    assert filtered.indexWithName("Na+(aq)")   < filtered.size()
    assert filtered.indexWithName("Cl-(aq)")   < filtered.size()
    assert filtered.indexWithName("NaCl(aq)")  < filtered.size()
    assert filtered.indexWithName("H2O(g)")    < filtered.size()

    filtered = specieslist.withElementsOf("HOCNaCl")

    assert filtered.size() == specieslist.size()
    assert filtered.indexWithName("H2O(aq)"  ) < filtered.size()
    assert filtered.indexWithName("H+(aq)"   ) < filtered.size()
    assert filtered.indexWithName("OH-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("H2(aq)"   ) < filtered.size()
    assert filtered.indexWithName("O2(aq)"   ) < filtered.size()
    assert filtered.indexWithName("Na+(aq)"  ) < filtered.size()
    assert filtered.indexWithName("Cl-(aq)"  ) < filtered.size()
    assert filtered.indexWithName("NaCl(aq)" ) < filtered.size()
    assert filtered.indexWithName("CO2(aq)"  ) < filtered.size()
    assert filtered.indexWithName("HCO3-(aq)") < filtered.size()
    assert filtered.indexWithName("CO3-2(aq)") < filtered.size()
    assert filtered.indexWithName("CH4(aq)"  ) < filtered.size()
    assert filtered.indexWithName("H2O(g)"   ) < filtered.size()
    assert filtered.indexWithName("CO2(g)"   ) < filtered.size()
    assert filtered.indexWithName("CH4(g)"   ) < filtered.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::append
    #-------------------------------------------------------------------------
    specieslist.append(Species("CaCO3(calcite)"))

    assert specieslist.indexWithName("CaCO3(calcite)") < specieslist.size()
    assert specieslist.indexWithFormula("CaCO3") < specieslist.size()

    #-------------------------------------------------------------------------
    # TESTING METHOD: SpeciesList::__iter__
    #-------------------------------------------------------------------------
    for i, species in enumerate(specieslist):
        assert species.name() == specieslist[i].name()
