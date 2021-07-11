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
    # Testing AqueousPhase with provided speciates and tags, so that species possessing them are excluded
    # -----------------------------------------------------------------------------------------------------------------

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C Na Cl"), exclude("organic")) )
    phases.add( GaseousPhase("CO2(g)") )

    phases = phases.convert()

    assert len(phases) == 2

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"

    phase0_species_list = ""
    for species in phases[0].species():
        phase0_species_list += species.name() + " "
    assert phase0_species_list[0:-1] == "H2O(aq) 2-Hydroxynonanoate- 2-Hydroxynonanoic(aq) CO(aq) CO2(aq) CO3-2 Cl- HClO(aq) ClO- ClO2- ClO3- ClO4- H+ H2(aq) HCO3- HO2- Nonanoate- Nonanoic-Acid(aq) Na+ NaCl(aq) NaOH(aq) O2(aq) OH- Nonanal(aq) H2O2(aq) HClO2(aq) HCl(aq)"

    phase1_species_list = ""
    for species in phases[0].species():
        phase1_species_list += species.name() + " "
    assert phase1_species_list[0:-1] == "CO2(g)"

    # AqueousPhase with provided tags, so that species possessing them are excluded

    phases = Phases(db)
    phases.add( AqueousPhase(exclude("organic")) )
    phases.add( GaseousPhase("CO2(g)") )

    phases = phases.convert()

    assert len(phases) == 2

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"

    phase0_species_list = ""
    for species in phases[0].species():
        phase0_species_list += species.name() + " "
    assert phase0_species_list[0:-1] == "H2O(aq) H+ H2(aq) HO2- O2(aq) OH- H2O2(aq)"

    phase1_species_list = ""
    for species in phases[0].species():
        phase1_species_list += species.name() + " "
    assert phase1_species_list[0:-1] == "CO2(g)"

    # Testing MineralPhases with provided tags, so that species possessing them are excluded

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C"), exclude("organic")) )
    phases.add( GaseousPhase("H2O(g) CO2(g)") )
    phases.add( MineralPhases(exclude("carbonate")) )

    phases = phases.convert()

    assert len(phases) == 3

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"
    assert phases[2].name() == "Graphite"

    phase0_species_list = ""
    for species in phases[0].species():
        phase0_species_list += species.name() + " "
    assert phase0_species_list[0:-1] == "H2O(aq) 2-Hydroxynonanoate- 2-Hydroxynonanoic(aq) CO(aq) CO2(aq) CO3-2 H+ H2(aq) HCO3- HO2- Nonanoate- Nonanoic-Acid(aq) O2(aq) OH- Nonanal(aq) H2O2(aq)"

    phase1_species_list = ""
    for species in phases[0].species():
        phase1_species_list += species.name() + " "
    assert phase1_species_list[0:-1] == "H2O(g) CO2(g)"

    phase2_species_list = ""
    for species in phases[0].species():
        phase2_species_list += species.name() + " "
    assert phase2_species_list[0:-1] == "Graphite"

    # Testing GaseousPhases with provided speciates and tags, so that species possessing them are excluded

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C"), exclude("organic")) )
    phases.add( GaseousPhase(speciate("H O C"), exclude("inert")) )
    phases.add( MineralPhases(exclude("carbonate")) )

    phases = phases.convert()

    assert len(phases) == 3

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"
    assert phases[2].name() == "Graphite"

    phase0_species_list = ""
    for species in phases[0].species():
        phase0_species_list += species.name() + " "
    assert phase0_species_list[0:-1] == "H2O(aq) 2-Hydroxynonanoate- 2-Hydroxynonanoic(aq) CO(aq) CO2(aq) CO3-2 H+ H2(aq) HCO3- HO2- Nonanoate- Nonanoic-Acid(aq) O2(aq) OH- Nonanal(aq) H2O2(aq)"

    phase1_species_list = ""
    for species in phases[0].species():
        phase1_species_list += species.name() + " "
    assert phase1_species_list[0:-1] == "CH4(g) C6H6O(g) o-Cresol(g) m-Cresol(g) p-Cresol(g) CO(g) CO2(g) C2H4(g) H2(g) H2O(g) O2(g)"

    phase2_species_list = ""
    for species in phases[0].species():
        phase2_species_list += species.name() + " "
    assert phase2_species_list[0:-1] == "Graphite"

    # Testing GaseousPhases with provided tags, so that species possessing them are excluded

    phases = Phases(db)
    phases.add( AqueousPhase(speciate("H O C Na Cl"), exclude("organic")) )
    phases.add( GaseousPhase(exclude("inert")) )
    phases.add( MineralPhases(exclude("carbonate")) )

    phases = phases.convert()

    assert len(phases) == 5

    assert phases[0].name() == "AqueousPhase"
    assert phases[1].name() == "GaseousPhase"
    assert phases[2].name() == "Graphite"
    assert phases[3].name() == "Halite"
    assert phases[4].name() == "Sodium-Oxide"

    phase0_species_list = ""
    for species in phases[0].species():
        phase0_species_list += species.name() + " "
    assert phase0_species_list[0:-1] == "H2O(aq) 2-Hydroxynonanoate- 2-Hydroxynonanoic(aq) CO(aq) CO2(aq) CO3-2 Cl- HClO(aq) ClO- ClO2- ClO3- ClO4- H+ H2(aq) HCO3- HO2- Nonanoate- Nonanoic-Acid(aq) Na+ NaCl(aq) NaOH(aq) O2(aq) OH- Nonanal(aq) H2O2(aq) HClO2(aq) HCl(aq)"

    phase1_species_list = ""
    for species in phases[0].species():
        phase1_species_list += species.name() + " "
    assert phase1_species_list[0:-1] == "CH4(g) C6H6O(g) o-Cresol(g) m-Cresol(g) p-Cresol(g) CO(g) CO2(g) C2H4(g) H2(g) H2O(g) O2(g)"

    phase2_species_list = ""
    for species in phases[0].species():
        phase2_species_list += species.name() + " "
    assert phase2_species_list[0:-1] == "Graphite"

    phase3_species_list = ""
    for species in phases[0].species():
        phase3_species_list += species.name() + " "
    assert phase3_species_list[0:-1] == "Halite"

    phase4_species_list = ""
    for species in phases[0].species():
        phase4_species_list += species.name() + " "
    assert phase4_species_list[0:-1] == "Sodium-Oxide"