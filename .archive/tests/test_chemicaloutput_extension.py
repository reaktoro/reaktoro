import pathlib
from typing import Tuple, Dict

import numpy as np
import pandas as pd
import pytest
from reaktoro import ChemicalEditor, ChemicalSystem, ReactionSystem, Partition, ChemicalState
from reaktoro import EquilibriumProblem, equilibrate, KineticPath


@pytest.fixture
def dict_with_properties_to_output():
    properties_and_components = {
        "time(units=minute)": "time(units=minute)",
        "pH": "pH",
        "Ca [mmolal]": "elementMolality(Ca units=mmolal)",
        "Mg [mmolal]": "elementMolality(Mg units=mmolal)",
        "Calcite [units=g]": "phaseMass(Calcite units=g)",
        "Dolomite [units=g]": "phaseMass(Dolomite units=g)"
    }
    return properties_and_components


@pytest.fixture
def brine_co2_path():
    editor = ChemicalEditor()
    editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3 MgCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")
    editor.addMineralPhase("Magnesite")
    editor.addMineralPhase("Dolomite")
    editor.addMineralPhase("Halite")

    editor.addMineralReaction("Calcite") \
        .setEquation("Calcite = Ca++ + CO3--") \
        .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, "cm2/g")

    editor.addMineralReaction("Magnesite") \
        .setEquation("Magnesite = Mg++ + CO3--") \
        .addMechanism("logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol") \
        .addMechanism("logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0") \
        .setSpecificSurfaceArea(10, "cm2/g")

    editor.addMineralReaction("Dolomite") \
        .setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--") \
        .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol") \
        .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5") \
        .setSpecificSurfaceArea(10, "cm2/g")

    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)

    partition = Partition(system)
    partition.setKineticSpecies(["Calcite", "Magnesite", "Dolomite"])

    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.setTemperature(60, "celsius")
    problem.setPressure(100, "bar")
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.5, "mol")
    problem.add("CO2", 1, "mol")

    state = equilibrate(problem)

    state.setSpeciesMass("Calcite", 100, "g")
    state.setSpeciesMass("Dolomite", 50, "g")

    path = KineticPath(reactions)
    path.setPartition(partition)

    return path, state


@pytest.fixture
def output_from_path(
    brine_co2_path: Tuple[KineticPath, ChemicalState],
    tmp_path: pathlib.Path,
    dict_with_properties_to_output: Dict,
):
    path, state = brine_co2_path
    output = path.output()
    output.filename(str(tmp_path / "test_output_path.txt"))
    for property_name, unit in dict_with_properties_to_output.items():
        if property_name == "time(units=minute)" or property_name == "pH":
            output.add(property_name)
        else:
            output.add(unit, property_name)

    t0, t1 = 0.0, 25.0
    path.solve(state, t0, t1, "hours")

    return output


def test_chemicaloutput_to_array(
        dict_with_properties_to_output, output_from_path
):
    output = output_from_path
    output_array = output.to_array()

    assert type(output_array) is np.ndarray
    assert output_array.shape[1] == len(dict_with_properties_to_output.keys())


def test_chemicaloutput_to_dict(
        dict_with_properties_to_output, output_from_path
):
    output = output_from_path
    output_dict = output.to_dict()

    assert type(output_dict) is dict
    assert output_dict.keys() == dict_with_properties_to_output.keys()


def test_chemicaloutput_to_dataframe(
        dict_with_properties_to_output, output_from_path
):
    output = output_from_path
    output_df = output.to_data_frame()

    assert type(output_df) is pd.DataFrame
    assert list(output_df.columns) == list(dict_with_properties_to_output.keys())
