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


def testAqueousProps():

    db = PhreeqcDatabase("phreeqc.dat")

    solution = AqueousPhase(speciate("H O Na Cl Ca C"))
    system = ChemicalSystem(db, solution)

    state = ChemicalState(system)
    state.temperature(25.0, "celsius")
    state.pressure(1.0, "bar")
    state.set("H2O" , 0.5000, "kg")
    state.set("Na+" , 0.1234, "mol")
    state.set("Cl-" , 0.2345, "mol")
    state.set("Ca+2", 0.3456, "mol")
    state.set("CO2" , 0.1234, "mol")

    aprops = AqueousProps(state)

    assert aprops.speciesMolality("Na+").val()  == pytest.approx(2.0*0.1234)
    assert aprops.speciesMolality("Cl-").val()  == pytest.approx(2.0*0.2345)
    assert aprops.speciesMolality("Ca+2").val() == pytest.approx(2.0*0.3456)
    assert aprops.speciesMolality("CO2").val()  == pytest.approx(2.0*0.1234)

    # TODO: Implement more tests here for AqueousProps -- see its C++ tests.


# The databases used for testing alkalinity calculations
tested_databases = [
    SupcrtDatabase("supcrt07"),
    SupcrtDatabase("supcrtbl"),
    PhreeqcDatabase("phreeqc.dat"),
    ThermoFunDatabase("aq17")
]

# The alkalinity factors in Wolf-Gladrow et al. (2007) used in method AqueousProps.alkalinity
alkfactors = {
    "Na+"   :  1.0,
    "Mg+2"  :  2.0,
    "Ca+2"  :  2.0,
    "K+"    :  1.0,
    "Sr+2"  :  2.0,
    "Cl-"   : -1.0,
    "Br-"   : -1.0,
    "NO3-"  : -1.0,
    "H3PO4" : -1.0,
    "H2PO4-": -1.0,
    "HPO4-2": -1.0,
    "PO4-3" : -1.0,
    "NH3"   :  1.0,
    "NH4+"  :  1.0,
    "SO4-2" : -2.0,
    "HSO4-" : -2.0,
    "F-"    : -1.0,
    "HF"    : -1.0,
    "NO2-"  : -1.0,
    "HNO2"  : -1.0,
}

# The aqueous compositions used for testing alkalinity calculations (species amounts may not be realistical on purpose; for testing purposes only)
tested_compositions = [
    [
        ("H2O"  , 1.0    , "kg") ,
        ("H+"   , 6.3e-9 , "mol"),
        ("OH-"  , 10e-6  , "mol"),
        ("Na+"  , 0.6    , "mol"),
        ("Cl-"  , 0.59767, "mol"),
        ("CO2"  , 8.3e-6 , "mol"),
        ("HCO3-", 1660e-6, "mol"),
        ("CO3-2", 331e-6 , "mol"),
    ],
    [
        ("H2O"  , 1.000, "kg") ,
        ("Na+"  , 1.234, "mol"),
        ("Cl-"  , 2.345, "mol"),
        ("CO2"  , 0.123, "mol"),
        ("HCO3-", 0.234, "mol"),
        ("CO3-2", 0.456, "mol"),
        ("Ca+2" , 0.578, "mol"),
        ("Mg+2" , 0.987, "mol"),
        ("NH3"  , 0.012, "mol"),
        ("NH4+" , 0.987, "mol"),
    ],
    [
        ("H2O"   , 0.5000, "kg") ,
        ("Na+"   , 0.0100, "mol"),
        ("Mg+2"  , 0.1200, "mol"),
        ("Ca+2"  , 0.2300, "mol"),
        ("K+"    , 0.3400, "mol"),
        ("Sr+2"  , 0.4500, "mol"),
        ("Cl-"   , 0.5600, "mol"),
        ("Br-"   , 0.6700, "mol"),
        ("NO3-"  , 0.7800, "mol"),
        ("H3PO4" , 0.8900, "mol"),
        ("H2PO4-", 0.9100, "mol"),
        ("HPO4-2", 0.1011, "mol"),
        ("PO4-3" , 0.1112, "mol"),
        ("NH3"   , 0.1213, "mol"),
        ("NH4+"  , 0.1314, "mol"),
        ("SO4-2" , 0.1415, "mol"),
        ("HSO4-" , 0.1516, "mol"),
        ("F-"    , 0.1617, "mol"),
        ("HF"    , 0.1718, "mol"),
        ("NO2-"  , 0.1819, "mol"),
        ("HNO2"  , 0.1920, "mol"),
    ],
]

@pytest.mark.parametrize("db", tested_databases)
@pytest.mark.parametrize("composition", tested_compositions)
def testAqueousPropsAlkalinity(db, composition):

    formulas = [entry[0] for entry in composition]

    solution = AqueousPhase(speciate(formulas))
    system = ChemicalSystem(db, solution)

    specieslist = system.species()

    state = ChemicalState(system)
    state.temperature(25.0, "celsius")
    state.pressure(1.0, "atm")

    for formula, quantity, unit in composition:
        idx = system.species().findWithFormula(formula)
        if idx < specieslist.size():
            state.set(idx, quantity, unit)

    props = ChemicalProps(state)
    aqprops = AqueousProps(state)

    expected_alk = 0.0
    for formula in formulas:
        idx = specieslist.findWithFormula(formula)
        if idx < specieslist.size():
            amount = state.speciesAmount(idx)
            alkfactor = alkfactors.get(formula, 0.0)
            expected_alk += alkfactor * amount

    liters = props.phaseProps(0).volume() * 1e3  # convert volume from m3 to L
    expected_alk /= liters

    assert aqprops.alkalinity().val() == pytest.approx(expected_alk)
