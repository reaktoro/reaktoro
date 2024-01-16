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
from pytest import approx

def testSurfaceAreaModelPower():

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

    model1 = SurfaceAreaModelPowerMolar("Calcite", 1000.0, 20.0, 3.0)
    model2 = SurfaceAreaModelPowerSpecific("Magnesite", 1000.0, 20.0, 3.0)
    model3 = SurfaceAreaModelPowerVolumetric("Quartz", 1000.0, 20.0, 3.0)

    assert model1(props).val() == approx(1.0)  # 1000 * (2 / 20)**3
    assert model2(props).val() == approx(1.0)  # 1000 * (2 / 20)**3
    assert model3(props).val() == approx(1.0)  # 1000 * (2 / 20)**3
