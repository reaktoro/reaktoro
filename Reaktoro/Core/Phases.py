# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2021 Allan Leal
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


def checkSpeicesInPhase(phase, expected_species):
    phase_species = ""
    for species in phase.species():
        phase_species += species.name() + " "
    assert phase_species[0:-1] == expected_species

def testPhases():

    db = Database()
    db.addSpecies( Species("H2O(aq)")                            )
    db.addSpecies( Species("H+")                                 )
    db.addSpecies( Species("OH-")                                )
    db.addSpecies( Species("H2(aq)")                             )
    db.addSpecies( Species("O2(aq)")                             )
    db.addSpecies( Species("Na+")                                )
    db.addSpecies( Species("Cl-")                                )
    db.addSpecies( Species("NaCl(aq)")                           )
    db.addSpecies( Species("HCl(aq)")                            )
    db.addSpecies( Species("NaOH(aq)")                           )
    db.addSpecies( Species("Ca++")                               )
    db.addSpecies( Species("Mg++")                               )
    db.addSpecies( Species("CO2(aq)")                            )
    db.addSpecies( Species("HCO3-")                              )
    db.addSpecies( Species("CO3--")                              )
    db.addSpecies( Species("CaCl2(aq)")                          )
    db.addSpecies( Species("MgCl2(aq)")                          )
    db.addSpecies( Species("SiO2(aq)")                           )
    db.addSpecies( Species("NaCl(s)").withName("Halite")         )
    db.addSpecies( Species("CaCO3(s)").withName("Calcite")       )
    db.addSpecies( Species("MgCO3(s)").withName("Magnesite")     )
    db.addSpecies( Species("CaMg(CO3)2(s)").withName("Dolomite") )
    db.addSpecies( Species("SiO2(s)").withName("Quartz")         )
    db.addSpecies( Species("CO2(g)")                             )
    db.addSpecies( Species("O2(g)")                              )
    db.addSpecies( Species("H2(g)")                              )
    db.addSpecies( Species("H2O(g)")                             )
    db.addSpecies( Species("CH4(g)")                             )
    db.addSpecies( Species("CO(g)")                              )

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O Na Cl C")) )
    phases.add( AqueousPhase("H2O(aq) H+ OH- Na+ Cl- HCO3- CO3-- CO2(aq)") )
    phases.add( GaseousPhase(speciate("H O C")) )
    phases.add( GaseousPhase("H2O(g) CO2(g)") )
    phases.add( MineralPhase("Halite") )
    phases.add( MineralPhases("Calcite Magnesite Dolomite Quartz") )

    phases = phases.convert()

    assert len(phases) == 9

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "AqueousPhase!"
    assert phases[2].name() == "GaseousPhase"
    assert phases[3].name() == "GaseousPhase!"
    assert phases[4].name() == "Halite"
    assert phases[5].name() == "Calcite"
    assert phases[6].name() == "Magnesite"
    assert phases[7].name() == "Dolomite"
    assert phases[8].name() == "Quartz"

    # -----------------------------------------------------------------------------------------------------------------
    # Testing Exclude tags functionality
    # -----------------------------------------------------------------------------------------------------------------

    db = Database()
    db.addSpecies( Species("H2O(aq)")                            )
    db.addSpecies( Species("H+")                                 )
    db.addSpecies( Species("OH-")                                )
    db.addSpecies( Species("H2(aq)")                             )
    db.addSpecies( Species("O2(aq)")                             )
    db.addSpecies( Species("Na+")                                )
    db.addSpecies( Species("Cl-")                                )
    db.addSpecies( Species("NaCl(aq)")                           )
    db.addSpecies( Species("HCl(aq)")                            )
    db.addSpecies( Species("NaOH(aq)")                           )
    db.addSpecies( Species("Ca++")                               )
    db.addSpecies( Species("Mg++")                               )
    db.addSpecies( Species("CO2(aq)")                            )
    db.addSpecies( Species("HCO3-")                              )
    db.addSpecies( Species("CO3--")                              )
    db.addSpecies( Species("CaCl2(aq)")                          )
    db.addSpecies( Species("MgCl2(aq)")                          )
    db.addSpecies( Species("SiO2(aq)")                           )
    db.addSpecies( Species("NaCl(s)").withName("Halite")         )
    db.addSpecies( Species("SiO2(s)").withName("Quartz")         )
    db.addSpecies( Species("CO2(g)")                             )
    db.addSpecies( Species("O2(g)")                              )
    db.addSpecies( Species("H2(g)")                              )
    db.addSpecies( Species("H2O(g)")                             )
    db.addSpecies( Species("CH4(g)")                             )
    db.addSpecies( Species("CO(g)")                              )

    db.addSpecies( Species("CaCO3(s)").withName("Calcite").withTags("carbonate")        )
    db.addSpecies( Species("MgCO3(s)").withName("Magnesite").withTags("carbonate")      )
    db.addSpecies( Species("CaMg(CO3)2(s)").withName("Dolomite").withTags("carbonate")  )
    db.addSpecies( Species("C(s)").withName("Graphite")                                 )
    db.addSpecies( Species("CaO(s)").withName("Lime")                                   )
    db.addSpecies( Species("N2(g)").withTags("inert")                                   )
    db.addSpecies( Species("C4H9OH").withTags("organic").withName("1-Butanol(aq)")      )
    db.addSpecies( Species("C4H8").withTags("organic").withName("1-Butene(aq)")         )
    db.addSpecies( Species("BaSO4(s)").withName("Barite").withTags("sulfate")           )
    db.addSpecies( Species("SrSO4(s)").withName("Celestite").withTags("sulfate")        )
    db.addSpecies( Species("PbSO4(s)").withName("Anglesite").withTags("sulfate")        )
    db.addSpecies( Species("CaSO4(s)").withName("Anhydrite").withTags("sulfate")        )

    # -----------------------------------------------------------------------------------------------------------------
    # Testing AqueousPhase with provided speciates
    # -----------------------------------------------------------------------------------------------------------------

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C Na Cl")) )
    phases.add( GaseousPhase("CO2(g)") )

    phases = phases.convert()

    assert len(phases) == 2

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"

    checkSpeicesInPhase(phases[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) CO2(aq) HCO3- CO3-- 1-Butanol(aq) 1-Butene(aq)")
    checkSpeicesInPhase(phases[1], "CO2(g)")

    # -----------------------------------------------------------------------------------------------------------------
    # Testing AqueousPhase with provided speciates and tags, so that species possessing them are excluded
    # -----------------------------------------------------------------------------------------------------------------

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C Na Cl"), exclude("organic")) )
    phases.add( GaseousPhase("CO2(g)") )

    phases = phases.convert()

    assert len(phases) == 2

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"

    checkSpeicesInPhase(phases[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) CO2(aq) HCO3- CO3--")
    checkSpeicesInPhase(phases[1], "CO2(g)")

    # -----------------------------------------------------------------------------------------------------------------
    # AqueousPhase with provided tags, so that species possessing them are excluded
    # -----------------------------------------------------------------------------------------------------------------

    phases = Phases(db)
    phases.add( AqueousPhase(exclude("organic")) )
    phases.add( GaseousPhase("CO2(g)") )

    phases = phases.convert()

    assert len(phases) == 2

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"

    checkSpeicesInPhase(phases[0], "H2O(aq) H+ OH- H2(aq) O2(aq)")
    checkSpeicesInPhase(phases[1], "CO2(g)")

    # -----------------------------------------------------------------------------------------------------------------
    # Testing MineralPhases with provided speciate symbols and tags, so that species possessing them are excluded
    # -----------------------------------------------------------------------------------------------------------------

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O"), exclude("organic")))
    phases.add( MineralPhases(speciate("C Ca O"), exclude("carbonate")) )

    phases = phases.convert()

    assert len(phases) == 3

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "Graphite"
    assert phases[2].name() == "Lime"

    checkSpeicesInPhase(phases[0], "H2O(aq) H+ OH- H2(aq) O2(aq)")
    checkSpeicesInPhase(phases[1], "Graphite")
    checkSpeicesInPhase(phases[2], "Lime")

    # -----------------------------------------------------------------------------------------------------------------
    # Testing MineralPhases with provided tags, so that species possessing them are excluded
    # -----------------------------------------------------------------------------------------------------------------

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C"), exclude("organic")) )
    phases.add( GaseousPhase("H2O(g) CO2(g)") )
    phases.add( MineralPhases(exclude("carbonate")) )

    phases = phases.convert()

    assert len(phases) == 3

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"
    assert phases[2].name() == "Graphite"

    checkSpeicesInPhase(phases[0], "H2O(aq) H+ OH- H2(aq) O2(aq) CO2(aq) HCO3- CO3--")
    checkSpeicesInPhase(phases[1], "H2O(g) CO2(g)")
    checkSpeicesInPhase(phases[2], "Graphite")

    # -----------------------------------------------------------------------------------------------------------------
    # Testing GaseousPhases with provided speciates and tags, so that species possessing them are excluded
    # -----------------------------------------------------------------------------------------------------------------

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C"), exclude("organic")) )
    phases.add( GaseousPhase(speciate("H O C N"), exclude("inert")) )
    phases.add( MineralPhases(exclude("carbonate")) )

    phases = phases.convert()

    assert len(phases) == 3

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"
    assert phases[2].name() == "Graphite"

    checkSpeicesInPhase(phases[0], "H2O(aq) H+ OH- H2(aq) O2(aq) CO2(aq) HCO3- CO3--")
    checkSpeicesInPhase(phases[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)")
    checkSpeicesInPhase(phases[2], "Graphite")

    # -----------------------------------------------------------------------------------------------------------------
    # Testing GaseousPhases with provided tags, so that species possessing them are excluded
    # -----------------------------------------------------------------------------------------------------------------

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C Na Cl"), exclude("organic")) )
    phases.add( GaseousPhase(exclude("inert")) )
    phases.add( MineralPhases(exclude("carbonate")) )

    phases = phases.convert()

    assert len(phases) == 4

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"
    assert phases[2].name() == "Halite"
    assert phases[3].name() == "Graphite"

    checkSpeicesInPhase(phases[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) CO2(aq) HCO3- CO3--")
    checkSpeicesInPhase(phases[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)")
    checkSpeicesInPhase(phases[2], "Halite")
    checkSpeicesInPhase(phases[3], "Graphite")
