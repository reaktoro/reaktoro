import numpy as np
import pandas as pd
import pytest

from collections import namedtuple
from python_tools import convert_dataframe_to_dict, convert_reaktoro_state_to_dict
from reaktoro import (
    ChemicalEditor,
    ChemicalSystem,
    Database,
    equilibrate,
    EquilibriumProblem,
    KineticPath,
    Partition,
    ReactionSystem,
)


@pytest.fixture(scope="session")
def kinect_problem_with_h2o_hcl_caco3_mgco3_co2_calcite():
    '''
    Build a kinetic problem with 1 kg of H2O, 1mmol of HCl which has calcite
    as a kinetic reaction
    '''

    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhase("H2O HCl CaCO3 MgCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")

    calciteReaction = editor.addMineralReaction("Calcite")
    calciteReaction.setEquation("Calcite = Ca++ + CO3--")
    calciteReaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
    calciteReaction.addMechanism(
        "logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0"
    )
    calciteReaction.setSpecificSurfaceArea(10, "cm2/g")

    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)

    partition = Partition(system)
    partition.setKineticPhases(["Calcite"])

    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.add("H2O", 1, "kg")
    problem.add("HCl", 1, "mmol")

    state = equilibrate(problem)

    state.setSpeciesMass("Calcite", 100, "g")

    return (state, reactions, partition)


@pytest.fixture(scope="session")
def kinect_problem_with_h2o_nacl_caco3_mgco3_hcl_co2_calcite_magnesite_dolomite_halite():

    '''
    Build a kinetic problem with 1 kg of H2O, 1 mol of NaCl and 1 mol of CO2
    which has the following kinetic reactions: calcite, Magnesite and Dolomite.
    '''
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)

    editor.addAqueousPhase("H2O NaCl CaCO3 MgCO3 HCl")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")
    editor.addMineralPhase("Magnesite")
    editor.addMineralPhase("Dolomite")
    editor.addMineralPhase("Halite")

    calciteReaction = editor.addMineralReaction("Calcite")
    calciteReaction.setEquation("Calcite = Ca++ + CO3--")
    calciteReaction.addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
    calciteReaction.addMechanism(
        "logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0"
    )
    calciteReaction.setSpecificSurfaceArea(10, "cm2/g")

    magnesiteReaction = editor.addMineralReaction("Magnesite")
    magnesiteReaction.setEquation("Magnesite = Mg++ + CO3--")
    magnesiteReaction.addMechanism("logk = -9.34 mol/(m2*s); Ea = 23.5 kJ/mol")
    magnesiteReaction.addMechanism(
        "logk = -6.38 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0"
    )
    magnesiteReaction.setSpecificSurfaceArea(10, "cm2/g")

    dolomiteReaction = editor.addMineralReaction("Dolomite")
    dolomiteReaction.setEquation("Dolomite = Ca++ + Mg++ + 2*CO3--")
    dolomiteReaction.addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol")
    dolomiteReaction.addMechanism(
        "logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5"
    )
    dolomiteReaction.setSpecificSurfaceArea(10, "cm2/g")

    system = ChemicalSystem(editor)
    reactions = ReactionSystem(editor)

    partition = Partition(system)
    partition.setKineticSpecies(["Calcite", "Magnesite", "Dolomite"])

    problem = EquilibriumProblem(system)
    problem.setPartition(partition)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 1, "mol")
    problem.add("CO2", 1, "mol")

    state = equilibrate(problem)

    state.setSpeciesMass("Calcite", 100, "g")
    state.setSpeciesMass("Dolomite", 50, "g")

    return (state, reactions, partition)


time_span = namedtuple("time_span", ["ti", "tf", "unit"])


@pytest.mark.serial
@pytest.mark.parametrize(
    "setup, time_span, checkedVariables",
    [
        (
            pytest.lazy_fixture(
                "kinect_problem_with_h2o_nacl_caco3_mgco3_hcl_co2_calcite_magnesite_dolomite_halite"
            ),
            time_span(0, 24, "hours"),
            [
                "time(units=hour)",
                "pH",
                "elementMolality(Ca units=molal)",
                "elementMolality(Mg units=molal)",
                "phaseMass(Calcite units=grams)",
                "phaseMass(Dolomite units=grams)",
            ],
        ),
        (
            pytest.lazy_fixture(
                "kinect_problem_with_h2o_nacl_caco3_mgco3_hcl_co2_calcite_magnesite_dolomite_halite"
            ),
            time_span(0, 48, "hours"),
            [
                "time(units=hour)",
                "pH",
                "elementMolality(Ca units=molal)",
                "elementMolality(Mg units=molal)",
                "phaseMass(Calcite units=grams)",
                "phaseMass(Dolomite units=grams)",
            ],
        ),
        (
            pytest.lazy_fixture(
                "kinect_problem_with_h2o_nacl_caco3_mgco3_hcl_co2_calcite_magnesite_dolomite_halite"
            ),
            time_span(0, 72, "hours"),
            [
                "time(units=hour)",
                "pH",
                "elementMolality(Ca units=molal)",
                "elementMolality(Mg units=molal)",
                "phaseMass(Calcite units=grams)",
                "phaseMass(Dolomite units=grams)",
            ],
        ),
        (
            pytest.lazy_fixture("kinect_problem_with_h2o_hcl_caco3_mgco3_co2_calcite"),
            time_span(0, 5, "minute"),
            [
                "time(units=minute)",
                "elementMolality(Ca units=mmolal)",
                "phaseMass(Calcite units=g)",
            ],
        ),
        (
            pytest.lazy_fixture("kinect_problem_with_h2o_hcl_caco3_mgco3_co2_calcite"),
            time_span(0, 10, "minute"),
            [
                "time(units=minute)",
                "elementMolality(Ca units=mmolal)",
                "phaseMass(Calcite units=g)",
            ],
        ),
        (
            pytest.lazy_fixture("kinect_problem_with_h2o_hcl_caco3_mgco3_co2_calcite"),
            time_span(0, 20, "minute"),
            [
                "time(units=minute)",
                "elementMolality(Ca units=mmolal)",
                "phaseMass(Calcite units=g)",
            ],
        ),
    ],
    ids=[
        "kinetic prob-h2o nacl caco3 mgco3 hcl co2 calcite magnesite dolomite halite 0 to 24 h",
        "kinetic prob-h2o nacl caco3 mgco3 hcl co2 calcite magnesite dolomite halite 0 to 48 h",
        "kinetic prob-h2o nacl caco3 mgco3 hcl co2 calcite magnesite dolomite halite 0 to 72 h",
        "kinetic prob-h2o hcl caco3 mgco3 co2 calcite 0 to 5 min",
        "kinetic prob-h2o hcl caco3 mgco3 co2 calcite 0 to 10 min",
        "kinetic prob-h2o hcl caco3 mgco3 co2 calcite 0 to 20 min",
    ],
)
def test_kinetic_path_solve_complete_path(
    num_regression, tmpdir, setup, time_span, checkedVariables
):
    """
    An integration test that checks result's reproducibility of 
    the calculation of a kinetic problem and check all the path 
    @param setup
        a tuple that has some objects from kineticProblemSetup.py
        (state, reactions, partition)
    @param time_span
        time information about the kinetic problem.
        time_span.ti = initial time
        time_span.tf = final time
        time_span.unit = ti and tf units
    @param checkedVariables 
        a list that has all the variables that will be tested
    """
    (state, reactions, partition) = setup

    path = KineticPath(reactions)

    path.setPartition(partition)

    output = path.output()
    output.filename(tmpdir.dirname + "/kinetictPathResult.txt")
    for checkedVariable in checkedVariables:
        output.add(checkedVariable)

    path.solve(state, time_span.ti, time_span.tf, time_span.unit)

    pathKineticTable = pd.read_csv(
        tmpdir.dirname + "/kinetictPathResult.txt",
        index_col=None,
        skiprows=1,
        delim_whitespace=True,
    )
    pathKineticTable.columns = checkedVariables

    pathKinectDic = convert_dataframe_to_dict(pathKineticTable)

    num_regression.check(pathKinectDic)


@pytest.mark.serial
@pytest.mark.parametrize(
    "setup, time_span",
    [
        (
            pytest.lazy_fixture(
                "kinect_problem_with_h2o_nacl_caco3_mgco3_hcl_co2_calcite_magnesite_dolomite_halite"
            ),
            time_span(0, 24, "hours"),
        ),
        (
            pytest.lazy_fixture(
                "kinect_problem_with_h2o_nacl_caco3_mgco3_hcl_co2_calcite_magnesite_dolomite_halite"
            ),
            time_span(0, 48, "hours"),
        ),
        (
            pytest.lazy_fixture(
                "kinect_problem_with_h2o_nacl_caco3_mgco3_hcl_co2_calcite_magnesite_dolomite_halite"
            ),
            time_span(0, 72, "hours"),
        ),
        (
            pytest.lazy_fixture("kinect_problem_with_h2o_hcl_caco3_mgco3_co2_calcite"),
            time_span(0, 5, "minute"),
        ),
        (
            pytest.lazy_fixture("kinect_problem_with_h2o_hcl_caco3_mgco3_co2_calcite"),
            time_span(0, 10, "minute"),
        ),
        (
            pytest.lazy_fixture("kinect_problem_with_h2o_hcl_caco3_mgco3_co2_calcite"),
            time_span(0, 20, "minute"),
        ),
    ],
    ids=[
        "kinetic prob-h2o nacl caco3 mgco3 hcl co2 calcite magnesite dolomite halite 0 to 24 h",
        "kinetic prob-h2o nacl caco3 mgco3 hcl co2 calcite magnesite dolomite halite 0 to 48 h",
        "kinetic prob-h2o nacl caco3 mgco3 hcl co2 calcite magnesite dolomite halite 0 to 72 h",
        "kinetic prob-h2o hcl caco3 mgco3 co2 calcite 0 to 5 min",
        "kinetic prob-h2o hcl caco3 mgco3 co2 calcite 0 to 10 min",
        "kinetic prob-h2o hcl caco3 mgco3 co2 calcite 0 to 20 min",
    ],
)
def test_kinetic_path_solve_final_state(num_regression, setup, time_span):
    """
    An integration test that checks result's reproducibility of 
    the calculation of a kinetic problem and only check the 
    final state 
    @param setup
        a tuple that has some objects from kineticProblemSetup.py
        (state, reactions, partition)
    @param time_span
        time information about the kinetic problem.
        time_span.ti = initial time
        time_span.tf = final time
        time_span.unit = ti and tf units
    """
    (state, reactions, partition) = setup

    path = KineticPath(reactions)

    path.setPartition(partition)

    path.solve(state, time_span.ti, time_span.tf, time_span.unit)

    stateDic = convert_reaktoro_state_to_dict(state)

    num_regression.check(stateDic)
