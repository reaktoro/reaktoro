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


def testMaterial():

    db = SupcrtDatabase("supcrtbl")

    solution = AqueousPhase(speciate("H O Na Cl Ca C Mg Si"), exclude("organic"))
    solution.setActivityModel(chain(ActivityModelHKF(), ActivityModelDrummond("CO2")))

    gases = GaseousPhase("H2O(g) CO2(g) O2(g) CH4(g)")
    gases.setActivityModel(ActivityModelPengRobinson())

    minerals = MineralPhases("Halite Calcite Magnesite Dolomite Quartz")

    system = ChemicalSystem(db, solution, gases, minerals)

    brine = Material(system)
    brine.add("H2O"    , 0.30, "kg")  # using formula
    brine.add("H2O(aq)", 0.70, "kg")  # using species name
    brine.add("NaCl"   , 1.00, "mol") # using formula
    brine.add("CaCl2"  , 0.10, "mol") # using formula
    brine.add("MgCl2"  , 0.05, "mol") # using formula

    rock = Material(system)
    rock.add("CaCO3", 1.00, "mol")    # using formula
    rock.add("SiO2" , 1.00, "mol")    # using formula

    mix = brine(1.0, "kg") + rock(1.0, "kg")

    bmix = mix.componentAmounts()  # the amounts of elements/charge in mix material

    state = ChemicalState(system)

    state = mix.equilibrate()

    assert mix.result().optima.succeeded
    assert mix.result().optima.iterations == 48

    assert state.temperature() == pytest.approx(25.0 + 273.15)
    assert state.pressure() == pytest.approx(1.0 * 1e5)
    assert state.componentAmounts().asarray() == pytest.approx(bmix)  # ensure computed state has element/charge amounts equal to those in mix material

    state = mix.equilibrate(60.0, "celsius", 10.0, "bar")

    assert mix.result().optima.succeeded
    assert mix.result().optima.iterations == 47

    assert state.temperature() == pytest.approx(60.0 + 273.15)
    assert state.pressure() == pytest.approx(10.0 * 1e5)
    assert state.componentAmounts().asarray() == pytest.approx(bmix)  # ensure computed state has element/charge amounts equal to those in mix material

