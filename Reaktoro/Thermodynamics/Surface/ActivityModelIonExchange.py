# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2021 Allan Leal
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
import numpy as np
import pytest


def initializeMoleFractions(species):

    idx = lambda formula: species.indexWithFormula(formula)

    n = 1e-6 * np.ones(species.size())

    n[idx("MgX2")] = 0.1
    n[idx("CaX2")] = 0.2
    n[idx("NaX")]  = 0.3

    return n / n.sum()


def testActivityModelIonExchage():

    db = PhreeqcDatabase("phreeqc.dat")

    # Expected species: AlOHX2 AlX3 BaX2 CaX2 CdX2 CuX2 FeX2 KX LiX MgX2 MnX2 NH4X NaX PbX2 SrX2 ZnX2
    species = db.species().withAggregateState(AggregateState.IonExchange).withCharge(0.0)

    # Check the names of the species
    assert species[0].name()  == "AlOHX2"
    assert species[1].name()  == "AlX3"
    assert species[2].name()  == "BaX2"
    assert species[3].name()  == "CaX2"
    assert species[4].name()  == "CdX2"
    assert species[5].name()  == "CuX2"
    assert species[6].name()  == "FeX2"
    assert species[7].name()  == "KX"
    assert species[8].name()  == "LiX"
    assert species[9].name()  == "MgX2"
    assert species[10].name() == "MnX2"
    assert species[11].name() == "NH4X"
    assert species[12].name() == "NaX"
    assert species[13].name() == "PbX2"
    assert species[14].name() == "SrX2"
    assert species[15].name() == "ZnX2"

    # Initialize input data for the ActivityProps
    T = autodiff.real(300.0)
    P = autodiff.real(12.3e5)
    x = initializeMoleFractions(species)

    # Construct the activity model function with the given ion exchange species.
    fn = ActivityModelIonExchangeGainesThomas()(species)

    # Create the ActivityProps object with the results.
    props = ActivityProps.create(species.size())

    # Evaluate the activity props function
    # fn(props, T, P, x)

    # Fetch exchanger equivalents in the ion exchange species
    #ze = props.extra.at(0)
    #
    # assert ze[0]  == 2
    # assert ze[1]  == 3
    # assert ze[2]  == 2
    # assert ze[3]  == 2
    # assert ze[4]  == 2
    # assert ze[5]  == 2
    # assert ze[6]  == 2
    # assert ze[7]  == 1
    # assert ze[8]  == 1
    # assert ze[9]  == 2
    # assert ze[10] == 2
    # assert ze[11] == 1
    # assert ze[12] == 1
    # assert ze[13] == 2
    # assert ze[14] == 2
    # assert ze[15] == 2

    # props.ln_a[0]  == Approx(-13.017)
    # props.ln_a[1]  == Approx(-12.6116)
    # props.ln_a[2]  == Approx(-13.017)
    # props.ln_a[3]  == Approx(-0.810957)
    # props.ln_a[4]  == Approx(-13.017)
    # props.ln_a[5]  == Approx(-13.017)
    # props.ln_a[6]  == Approx(-13.017)
    # props.ln_a[7]  == Approx(-13.7102)
    # props.ln_a[8]  == Approx(-13.7102)
    # props.ln_a[9]  == Approx(-1.5041)
    # props.ln_a[10] == Approx(-13.017)
    # props.ln_a[11] == Approx(-13.7102)
    # props.ln_a[12] == Approx(-1.09864)
    # props.ln_a[13] == Approx(-13.017)
    # props.ln_a[14] == Approx(-13.017)
    # props.ln_a[15] == Approx(-13.017)
