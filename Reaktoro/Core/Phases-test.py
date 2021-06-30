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

