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


def testThermoFunDatabase():
    aq17 = ThermoFunDatabase("aq17")

    T = 298.15
    P = 1.0e5

    #-------------------------------------------------------------------
    # Testing attributes and thermodynamic properties of H2O@
    #-------------------------------------------------------------------
    species = aq17.species().get("H2O@")
    assert species.formula().equivalent("H2O")
    assert species.substance() == "Water HGK"
    assert species.aggregateState() == AggregateState.Aqueous
    assert species.charge() == 0
    assert species.molarMass() == pytest.approx(0.0180153)

    props = species.standardThermoProps(T, P)
    assert props.G0[0]  == pytest.approx(-2.371817e+05)
    assert props.H0[0]  == pytest.approx(-2.858310e+05)
    assert props.V0[0]  == pytest.approx( 1.806862e-05)
    assert props.Cp0[0] == pytest.approx( 7.532758e+01)

    #-------------------------------------------------------------------
    # Testing attributes and thermodynamic properties of CO3-2
    #-------------------------------------------------------------------
    species = aq17.species().get("CO3-2")
    assert species.formula().equivalent("CO3-2")
    assert species.substance() == "CO3-2 carbonate ion"
    assert species.aggregateState() == AggregateState.Aqueous
    assert species.charge() == -2
    assert species.molarMass() == pytest.approx(0.0600100979)

    props = species.standardThermoProps(T, P)
    assert props.G0[0]  == pytest.approx(-5.279830e+05)
    assert props.H0[0]  == pytest.approx(-6.752359e+05)
    assert props.V0[0]  == pytest.approx(-6.063738e-06)
    assert props.Cp0[0] == pytest.approx(-3.228612e+02)

    #-------------------------------------------------------------------
    # Testing attributes and thermodynamic properties of Ca+2
    #-------------------------------------------------------------------
    species = aq17.species().get("Ca+2")
    assert species.formula().equivalent("Ca+2")
    assert species.substance() == "Ca+2 ion"
    assert species.aggregateState() == AggregateState.Aqueous
    assert species.charge() == +2
    assert species.molarMass() == pytest.approx(0.040076902)

    props = species.standardThermoProps(T, P)
    assert props.G0[0]  == pytest.approx(-5.528210e+05)
    assert props.H0[0]  == pytest.approx(-5.431003e+05)
    assert props.V0[0]  == pytest.approx(-1.844093e-05)
    assert props.Cp0[0] == pytest.approx(-3.099935e+01)

    #-------------------------------------------------------------------
    # Testing attributes and thermodynamic properties of CO2
    #-------------------------------------------------------------------
    species = aq17.species().get("CO2")
    assert species.formula().equivalent("CO2")
    assert species.substance() == "Carbon dioxide (CO2)"
    assert species.aggregateState() == AggregateState.Gas
    assert species.charge() == 0
    assert species.molarMass() == pytest.approx(0.0440096006)

    props = species.standardThermoProps(T, P)
    assert props.G0[0]  == pytest.approx(-3.943510e+05)
    assert props.H0[0]  == pytest.approx(-3.935472e+05)
    assert props.V0[0]  == pytest.approx( 2.466434e-02)
    assert props.Cp0[0] == pytest.approx( 3.710812e+01)

    #-------------------------------------------------------------------
    # Testing attributes and thermodynamic properties of Calcite
    #-------------------------------------------------------------------
    species = aq17.species().get("Calcite")
    assert species.formula().equivalent("CaCO3")
    assert species.substance() == "Calcite (cc)"
    assert species.aggregateState() == AggregateState.CrystallineSolid
    assert species.charge() == 0
    assert species.molarMass() == pytest.approx(0.1000869999)

    props = species.standardThermoProps(T, P)
    assert props.G0[0]  == pytest.approx(-1.129195e+06)
    assert props.H0[0]  == pytest.approx(-1.207470e+06)
    assert props.V0[0]  == pytest.approx( 3.689000e-05)
    assert props.Cp0[0] == pytest.approx( 8.337073e+01)

    with pytest.raises(RuntimeError):
        assert ThermoFunDatabase("not-a-valid-file-name")

    with pytest.raises(RuntimeError):
        assert ThermoFunDatabase.withName("not-a-valid-file-name")

    with pytest.raises(RuntimeError):
        assert ThermoFunDatabase.fromFile("not-a-valid-file-name")
