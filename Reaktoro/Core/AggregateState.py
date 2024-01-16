# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
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


def testAggregateState():
    assert parseAggregateState("g")   == AggregateState.Gas
    assert parseAggregateState("l")   == AggregateState.Liquid
    assert parseAggregateState("s")   == AggregateState.Solid
    assert parseAggregateState("pl")  == AggregateState.Plasma
    assert parseAggregateState("cd")  == AggregateState.CondensedPhase
    assert parseAggregateState("fl")  == AggregateState.Fluid
    assert parseAggregateState("lc")  == AggregateState.LiquidCrystal
    assert parseAggregateState("cr")  == AggregateState.CrystallineSolid
    assert parseAggregateState("am")  == AggregateState.AmorphousSolid
    assert parseAggregateState("vit") == AggregateState.Vitreous
    assert parseAggregateState("ads") == AggregateState.Adsorbed
    assert parseAggregateState("mon") == AggregateState.Monomeric
    assert parseAggregateState("pol") == AggregateState.Polymeric
    assert parseAggregateState("ss")  == AggregateState.SolidSolution
    assert parseAggregateState("ex")  == AggregateState.IonExchange
    assert parseAggregateState("aq")  == AggregateState.Aqueous
    assert parseAggregateState("xy")  == AggregateState.Undefined

    assert parseAggregateState("Gas")              == AggregateState.Gas
    assert parseAggregateState("Liquid")           == AggregateState.Liquid
    assert parseAggregateState("Solid")            == AggregateState.Solid
    assert parseAggregateState("Plasma")           == AggregateState.Plasma
    assert parseAggregateState("CondensedPhase")   == AggregateState.CondensedPhase
    assert parseAggregateState("Fluid")            == AggregateState.Fluid
    assert parseAggregateState("LiquidCrystal")    == AggregateState.LiquidCrystal
    assert parseAggregateState("CrystallineSolid") == AggregateState.CrystallineSolid
    assert parseAggregateState("AmorphousSolid")   == AggregateState.AmorphousSolid
    assert parseAggregateState("Vitreous")         == AggregateState.Vitreous
    assert parseAggregateState("Adsorbed")         == AggregateState.Adsorbed
    assert parseAggregateState("Monomeric")        == AggregateState.Monomeric
    assert parseAggregateState("Polymeric")        == AggregateState.Polymeric
    assert parseAggregateState("SolidSolution")    == AggregateState.SolidSolution
    assert parseAggregateState("IonExchange")      == AggregateState.IonExchange
    assert parseAggregateState("Aqueous")          == AggregateState.Aqueous
    assert parseAggregateState("Undefined")        == AggregateState.Undefined
    assert parseAggregateState("XYZ")              == AggregateState.Undefined

    assert identifyAggregateState("XYZ(g)")       == AggregateState.Gas
    assert identifyAggregateState("XYZ(l)")       == AggregateState.Liquid
    assert identifyAggregateState("XYZ(s)")       == AggregateState.Solid
    assert identifyAggregateState("XYZ(s, xyz)")  == AggregateState.Solid
    assert identifyAggregateState("XYZ(pl)")      == AggregateState.Plasma
    assert identifyAggregateState("XYZ(cd)")      == AggregateState.CondensedPhase
    assert identifyAggregateState("XYZ(fl)")      == AggregateState.Fluid
    assert identifyAggregateState("XYZ(lc)")      == AggregateState.LiquidCrystal
    assert identifyAggregateState("XYZ(cr)")      == AggregateState.CrystallineSolid
    assert identifyAggregateState("XYZ(am)")      == AggregateState.AmorphousSolid
    assert identifyAggregateState("XYZ(vit)")     == AggregateState.Vitreous
    assert identifyAggregateState("XYZ(ads)")     == AggregateState.Adsorbed
    assert identifyAggregateState("XYZ(mon)")     == AggregateState.Monomeric
    assert identifyAggregateState("XYZ(pol)")     == AggregateState.Polymeric
    assert identifyAggregateState("XYZ(ss)")      == AggregateState.SolidSolution
    assert identifyAggregateState("XYZ(ex)")      == AggregateState.IonExchange
    assert identifyAggregateState("XYZ(aq)")      == AggregateState.Aqueous
    assert identifyAggregateState("XYZ(aq, uvw)") == AggregateState.Aqueous

    assert identifyAggregateState("XYZ-")         == AggregateState.Aqueous
    assert identifyAggregateState("XYZ--")        == AggregateState.Aqueous
    assert identifyAggregateState("XYZ---")       == AggregateState.Aqueous
    assert identifyAggregateState("XYZ-2")        == AggregateState.Aqueous
    assert identifyAggregateState("XYZ-3")        == AggregateState.Aqueous

    assert identifyAggregateState("XYZ+")         == AggregateState.Aqueous
    assert identifyAggregateState("XYZ++")        == AggregateState.Aqueous
    assert identifyAggregateState("XYZ+++")       == AggregateState.Aqueous
    assert identifyAggregateState("XYZ+2")        == AggregateState.Aqueous
    assert identifyAggregateState("XYZ+3")        == AggregateState.Aqueous
    assert identifyAggregateState("XYZ[2-]")      == AggregateState.Aqueous

    assert identifyAggregateState("XYZ-(pl)")     == AggregateState.Plasma
    assert identifyAggregateState("XYZ--(pl)")    == AggregateState.Plasma
    assert identifyAggregateState("XYZ---(pl)")   == AggregateState.Plasma
    assert identifyAggregateState("XYZ-2(pl)")    == AggregateState.Plasma
    assert identifyAggregateState("XYZ-3(pl)")    == AggregateState.Plasma

    assert identifyAggregateState("XYZ+(pl)")     == AggregateState.Plasma
    assert identifyAggregateState("XYZ++(pl)")    == AggregateState.Plasma
    assert identifyAggregateState("XYZ+++(pl)")   == AggregateState.Plasma
    assert identifyAggregateState("XYZ[3+](pl)")  == AggregateState.Plasma
    assert identifyAggregateState("XYZ+2(pl)")    == AggregateState.Plasma
    assert identifyAggregateState("XYZ+3(pl)")    == AggregateState.Plasma

    assert identifyAggregateState("XYZ")          == AggregateState.Undefined
    assert identifyAggregateState("XYZ(xy)")      == AggregateState.Undefined
