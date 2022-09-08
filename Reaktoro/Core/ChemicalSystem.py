# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2022 Allan Leal
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


def testChemicalSystem():
    db = SupcrtDatabase("supcrtbl")

    phases = Phases(db)
    phases.add( AqueousPhase("H2O(aq) H+ OH- Na+ Cl- Ca+2 Mg+2 HCO3- CO3-2 CO2(aq) SiO2(aq)") )
    phases.add( GaseousPhase("H2O(g) CO2(g)") )
    phases.add( MineralPhase("Halite") )
    phases.add( MineralPhase("Calcite") )
    phases.add( MineralPhase("Magnesite") )
    phases.add( MineralPhase("Dolomite") )
    phases.add( MineralPhase("Quartz") )

    system = ChemicalSystem(phases)

    assert ChemicalSystem().id() == system.id() + 1  # new ChemicalSystem object has id = id0 + 1 and original system continues to have id = id0
    assert ChemicalSystem().id() == system.id() + 2  # new ChemicalSystem object has id = id0 + 2 and original system continues to have id = id0

    assert system.database().elements().size() == db.elements().size()
    assert system.database().species().size() == db.species().size()

    assert system.elements().size() == 8
    assert system.species().size() == 18
    assert system.phases().size() == 7

    assert system.element(0).symbol() == "H"
    assert system.element(1).symbol() == "C"

    assert system.species(0).name() == "H2O(aq)"
    assert system.species(1).name() == "H+"

    assert system.phase(0).name() == "AqueousPhase"
    assert system.phase(1).name() == "GaseousPhase"

    assert (system.formulaMatrix() ==
        [[ 2,  1,  1,  0,  0,  0,  0,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0],
         [ 0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  1,  0,  1,  1,  2,  0],
         [ 1,  0,  1,  0,  0,  0,  0,  3,  3,  2,  2,  1,  2,  0,  3,  3,  6,  2],
         [ 0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0],
         [ 0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0],
         [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1],
         [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0],
         [ 0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0],
         [ 0,  1, -1,  1, -1,  2,  2, -1, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0]]).all()


    system = ChemicalSystem(db,
        AqueousPhase("H2O(aq) H+ OH- Na+ Cl- Ca+2 Mg+2 HCO3- CO3-2 CO2(aq) SiO2(aq)"),
        GaseousPhase("H2O(g) CO2(g)"),
        MineralPhase("Halite"),
        MineralPhase("Calcite"),
        MineralPhase("Magnesite"),
        MineralPhase("Dolomite"),
        MineralPhase("Quartz")
    )

    assert system.elements().size() == 8
    assert system.species().size() == 18
    assert system.phases().size() == 7
