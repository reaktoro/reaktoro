# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
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

def test_ReactionParams():
    reaction = ReactionParams()
    analytic = reaction.analytic
    assert not analytic
    analytic.append(1)
    assert analytic
    assert reaction.analytic
    assert reaction.analytic[0] == 1
    assert analytic is reaction.analytic


def test_SpeciesThermoData():
    species = SpeciesThermoData()
    assert species.properties is None
    assert species.reaction is None
    assert species.phreeqc is None

    species.properties = SpeciesThermoInterpolatedProperties()
    species.reaction = ReactionThermoInterpolatedProperties()
    species.phreeqc = SpeciesThermoParamsPhreeqc()
    species.phreeqc.reaction = ReactionParams()
    assert species.properties is not None
    assert species.reaction is not None
    assert species.phreeqc is not None
    assert species.phreeqc.reaction is not None


def test_AqueousSpeciesThermoData():
    aqueous_species = AqueousSpeciesThermoData()
    assert aqueous_species.hkf is None
    aqueous_species.hkf = AqueousSpeciesThermoParamsHKF()
    assert aqueous_species.hkf is not None


def test_GaseousSpeciesThermoData():
    gaseous_species = GaseousSpeciesThermoData()
    assert gaseous_species.hkf is None
    gaseous_species.hkf = GaseousSpeciesThermoParamsHKF()
    assert gaseous_species.hkf is not None


def test_MineralSpeciesThermoData():
    mineral_species = MineralSpeciesThermoData()
    assert mineral_species.hkf is None
    mineral_species.hkf = MineralSpeciesThermoParamsHKF()
    assert mineral_species.hkf is not None

