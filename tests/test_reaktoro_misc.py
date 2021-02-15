import pytest
from reaktoro import *


@pytest.fixture()
def create_minimum_system():
    """
    Build a system with only a aqueous phase
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H2O(l)")

    system = ChemicalSystem(editor)

    return system


def check_dict(chemical):
    dict1 = {'test':chemical}
    dict1[chemical] = 'test'

    dict2 = {chemical:'test'}
    dict2['test'] = chemical
    assert dict1['test'] == dict2['test'] and dict1[chemical] == dict2[chemical]


def test_dict_on_chemical_state(create_minimum_system):
    state = ChemicalState(create_minimum_system)
    check_dict(state)


def test_dict_on_chemical_quantity(create_minimum_system):
    quantity = ChemicalQuantity(create_minimum_system)
    check_dict(quantity)


def test_VectorDouble():
    vector = VectorDouble()
    assert not vector
    vector.append(15)
    vector.append(10)
    assert vector[0] == 15
    assert vector[1] == 10
    vector[0] = 1
    assert vector[0] == 1
    vector.pop()
    assert vector
    vector.pop()
    assert not vector
    vector.extend(VectorDouble(range(10)))
    assert len(vector) == 10
    for i in range(10):
        assert vector[i] == i


def test_ThermoScalar():
    T_input = 273.15
    T = ThermoScalar(T_input)

    assert T.val == T_input
    assert T.ddT == 0.0
    assert T.ddP == 0.0

    T_ddT = 10.0
    T_ddP = 30.0
    T_full = ThermoScalar(T_input, T_ddT, T_ddP)

    assert T_full.val == T_input
    assert T_full.ddT == T_ddT
    assert T_full.ddP == T_ddP
