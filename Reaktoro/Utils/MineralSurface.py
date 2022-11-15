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
from pytest import approx


def testMineralSurface():

    db = SupcrtDatabase("supcrtbl")

    system = ChemicalSystem(db,
        AqueousPhase("H2O(aq) H+ OH- Ca+2 HCO3- CO2(aq) CO3-2"),
        MineralPhase("Calcite"),
        MineralPhase("Magnesite"),
        MineralPhase("Quartz")
    )

    state = ChemicalState(system)

    state.set("H2O(aq)", 1.0, "kg")
    state.scalePhaseAmount("Calcite", 2.0, "mol")
    state.scalePhaseMass("Magnesite", 2.0, "kg")
    state.scalePhaseVolume("Quartz" , 2.0, "m3")

    props = ChemicalProps(state)

    # When MineralSurface is set with a constant surface area model
    surface1 = MineralSurface("Calcite", 1.23, "m2")
    surface2 = MineralSurface("Magnesite", 1.23, "cm2")
    surface3 = MineralSurface("Quartz", 1.23, "mm2")

    assert surface1.name() == "Calcite"
    assert surface2.name() == "Magnesite"
    assert surface3.name() == "Quartz"

    assert surface1.areaModel()(props).val() == approx(1.23)
    assert surface2.areaModel()(props).val() == approx(1.23e-4)
    assert surface3.areaModel()(props).val() == approx(1.23e-6)

    # When MineralSurface is set with a linear surface area model
    surface1 = MineralSurface("Calcite", 0.5, "cm2/mmol")
    surface2 = MineralSurface("Magnesite", 0.5, "cm2/g")
    surface3 = MineralSurface("Quartz", 0.5, "mm2/mm3")

    assert surface1.name() == "Calcite"
    assert surface2.name() == "Magnesite"
    assert surface3.name() == "Quartz"

    assert surface1.areaModel()(props).val() == approx(1.0e-1)
    assert surface2.areaModel()(props).val() == approx(1.0e-1)
    assert surface3.areaModel()(props).val() == approx(1.0e+3)

    # When MineralSurface is set with a power law surface area model
    surface1 = MineralSurface("Calcite", 1000.0, "m2", 20.0e6, "umol", 3.0)
    surface2 = MineralSurface("Magnesite", 1000.0e4, "cm2", 20.0e3, "g", 3.0)
    surface3 = MineralSurface("Quartz", 1000.0e6, "mm2", 20.0e9, "mm3", 3.0)

    assert surface1.name() == "Calcite"
    assert surface2.name() == "Magnesite"
    assert surface3.name() == "Quartz"

    assert surface1.areaModel()(props).val() == approx(1.0)  # 1000 * (2 / 20)**3
    assert surface2.areaModel()(props).val() == approx(1.0)  # 1000 * (2 / 20)**3
    assert surface3.areaModel()(props).val() == approx(1.0)  # 1000 * (2 / 20)**3

