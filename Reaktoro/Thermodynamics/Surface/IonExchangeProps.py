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
import pytest
import numpy as np

def testIonExchangeProps():

    db = PhreeqcDatabase("phreeqc.dat")

    # Define an aqueous phase
    solution = AqueousPhase(speciate("H O C Ca Na Mg Cl"))
    solution.setActivityModel(ActivityModelHKF())

    # Define an ion exchange phase
    exchange = IonExchangePhase("NaX CaX2 KX AlX3 MgX2")
    exchange.setActivityModel(ActivityModelIonExchange())

    # Create chemical system
    system = ChemicalSystem(db, solution, exchange)

    T = 25.0 # temperature in celsius
    P = 1.0  # pressure in bar

    # Define equilibrium solver
    solver = EquilibriumSolver(system)

    # Initialize ion exchange props
    exprops = IonExchangeProps(system)

    # Fetch the phase from the ion exchange props
    phase = exprops.phase()

    # Auxiliary variables
    species = system.species()
    exspecies = phase.species()
    exelements = phase.elements()
    num_species = species.size()

    # ------------------------------------------------------------------------
    # Testing correct initialization of the `IonExchangeProps` instance
    # ------------------------------------------------------------------------

    assert phase.species().size()   == 5
    assert phase.elements().size()  == 6
    assert phase.species(0).name()  == "NaX"
    assert phase.species(1).name()  == "CaX2"
    assert phase.species(2).name()  == "KX"
    assert phase.species(3).name()  == "AlX3"
    assert phase.species(4).name()  == "MgX2"

    # ------------------------------------------------------------------------
    # Testing when species have zero amounts
    # ------------------------------------------------------------------------

    n = np.zeros(num_species)
    state = ChemicalState(system)
    state.setTemperature(T, "celsius")
    state.setPressure(P, "bar")
    state.setSpeciesAmounts(n)

    with pytest.raises(Exception):
        exprops.update(state)

    # ------------------------------------------------------------------------
    #  Testing when species have nonzero amounts
    # ------------------------------------------------------------------------

    n = np.ones(num_species)
    state = ChemicalState(system)
    state.setTemperature(T, "celsius")
    state.setPressure(P, "bar")
    state.setSpeciesAmounts(n)

    # Update ion exchange properties
    exprops.update(state)

    # Check amounts of all species
    for i in range(0, exspecies.size()):
        assert exprops.speciesAmounts()[i] == pytest.approx(1.00)

    # Check equivalences of all species
    assert exprops.speciesEquivalence("NaX"  )[0] == pytest.approx(1.0)
    assert exprops.speciesEquivalence("CaX2" )[0] == pytest.approx(2.0)
    assert exprops.speciesEquivalence("KX"   )[0] == pytest.approx(1.0)
    assert exprops.speciesEquivalence("AlX3" )[0] == pytest.approx(3.0)
    assert exprops.speciesEquivalence("MgX2" )[0] == pytest.approx(2.0)

    # Check equivalent fractions of all species
    assert exprops.speciesEquivalentFraction("NaX"  )[0] == pytest.approx(0.1111111111111111)
    assert exprops.speciesEquivalentFraction("CaX2" )[0] == pytest.approx(0.2222222222222222)
    assert exprops.speciesEquivalentFraction("KX"   )[0] == pytest.approx(0.1111111111111111)
    assert exprops.speciesEquivalentFraction("AlX3" )[0] == pytest.approx(0.3333333333333333)
    assert exprops.speciesEquivalentFraction("MgX2" )[0] == pytest.approx(0.2222222222222222)

    # Check amounts of all elements
    assert exprops.elementAmount("Al")[0] == pytest.approx(1.0)
    assert exprops.elementAmount("Ca")[0] == pytest.approx(1.0)
    assert exprops.elementAmount("K" )[0] == pytest.approx(1.0)
    assert exprops.elementAmount("Na")[0] == pytest.approx(1.0)
    assert exprops.elementAmount("Mg")[0] == pytest.approx(1.0)
    assert exprops.elementAmount("X" )[0] == pytest.approx(9.0)

    # ------------------------------------------------------------------------
    # Testing when state is a brine"
    # ------------------------------------------------------------------------

    # Define initial equilibrium state
    state = ChemicalState(system)
    state.setTemperature(T, "celsius")
    state.setPressure(P, "bar")
    state.setSpeciesMass("H2O"   , 1.00, "kg")
    state.setSpeciesAmount("Na+" , 1.00, "mmol")
    state.setSpeciesAmount("Ca+2", 1.00, "mmol")
    state.setSpeciesAmount("Mg+2", 1.00, "mmol")
    state.setSpeciesAmount("NaX" , 1.00, "umol")
    solver.solve(state)

    # Update ion exchange properties
    exprops.update(state)

    # Check amounts of all species
    assert exprops.speciesAmount("NaX" )[0] == pytest.approx(9.84067523830254e-09  )
    assert exprops.speciesAmount("CaX2")[0] == pytest.approx(3.03996555695886e-07  )
    assert exprops.speciesAmount("KX"  )[0] == pytest.approx(1e-16                 )
    assert exprops.speciesAmount("AlX3")[0] == pytest.approx(1e-16                 )
    assert exprops.speciesAmount("MgX2")[0] == pytest.approx(1.9108310668496267e-07)

    # Check equivalences of all species
    assert exprops.speciesEquivalence("NaX"  )[0] == pytest.approx(9.84067523830254e-09  )
    assert exprops.speciesEquivalence("CaX2" )[0] == pytest.approx(6.07993111391772e-07  )
    assert exprops.speciesEquivalence("KX"   )[0] == pytest.approx(1e-16       )
    assert exprops.speciesEquivalence("AlX3" )[0] == pytest.approx(3e-16       )
    assert exprops.speciesEquivalence("MgX2" )[0] == pytest.approx(3.8216621336992535e-07)

    # Check equivalent fractions of all species
    assert exprops.speciesEquivalentFraction("NaX"  )[0] == pytest.approx(0.00984067523436627)
    assert exprops.speciesEquivalentFraction("CaX2" )[0] == pytest.approx(0.6079931111485748 )
    assert exprops.speciesEquivalentFraction("KX"   )[0] == pytest.approx(1e-10     )
    assert exprops.speciesEquivalentFraction("AlX3" )[0] == pytest.approx(3e-10     )
    assert exprops.speciesEquivalentFraction("MgX2" )[0] == pytest.approx(0.38216621321705885)

    # Check log10 gammas of all species
    assert exprops.speciesActivityCoefficientLg("NaX"  )[0] == pytest.approx(-0.0311024451808121 )
    assert exprops.speciesActivityCoefficientLg("CaX2" )[0] == pytest.approx(-0.12284266932733671)
    assert exprops.speciesActivityCoefficientLg("KX"   )[0] == pytest.approx(-0.03177765462046158)
    assert exprops.speciesActivityCoefficientLg("AlX3" )[0] == pytest.approx(-0.2575985613718797 )
    assert exprops.speciesActivityCoefficientLg("MgX2" )[0] == pytest.approx(-0.12146978806964734)

    # Check amounts of all elements
    assert exprops.elementAmount("Na")[0] == pytest.approx(9.84067523830254e-09  )
    assert exprops.elementAmount("Mg")[0] == pytest.approx(1.9108310668496267e-07)
    assert exprops.elementAmount("Al")[0] == pytest.approx(1e-16                 )
    assert exprops.elementAmount("K" )[0] == pytest.approx(1e-16                 )
    assert exprops.elementAmount("Ca")[0] == pytest.approx(3.03996555695886e-07  )
    assert exprops.elementAmount("X" )[0] == pytest.approx(1.0000000004e-06      )

    # Test convenience methods species and element amounts
    for s in exspecies:
        name = s.name()
        idx = exspecies.index(name)

        assert exprops.speciesAmount(name)[0]             == pytest.approx(exprops.speciesAmounts()[idx])
        assert exprops.speciesEquivalence(name)[0]        == pytest.approx(exprops.speciesEquivalences()[idx])
        assert exprops.speciesEquivalentFraction(name)[0] == pytest.approx(exprops.speciesEquivalentFractions()[idx])
        assert exprops.speciesActivityCoefficientLg(name)[0]         == pytest.approx(exprops.speciesActivityCoefficientsLg()[idx])

        assert exprops.speciesAmount(idx)[0]             == pytest.approx(exprops.speciesAmounts()[idx])
        assert exprops.speciesEquivalence(idx)[0]        == pytest.approx(exprops.speciesEquivalences()[idx])
        assert exprops.speciesEquivalentFraction(idx)[0] == pytest.approx(exprops.speciesEquivalentFractions()[idx])
        assert exprops.speciesActivityCoefficientLg(idx)[0]         == pytest.approx(exprops.speciesActivityCoefficientsLg()[idx])



    for e in exelements:
        symbol = e.symbol()
        idx = exelements.index(symbol)

        assert exprops.elementAmount(symbol)[0] == pytest.approx(exprops.elementAmounts()[idx])
        assert exprops.elementAmount(idx)[0] == pytest.approx(exprops.elementAmounts()[idx])
