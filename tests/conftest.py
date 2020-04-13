# -*- coding: utf-8 -*-

import pytest
import numpy as np

from reaktoro import ChemicalEditor, ChemicalSystem, Partition

from python_tools import convert_table_to_dict, convert_reaktoro_state_to_dict
from reaktoro import (
    Database,
    ChemicalEditor,
    ChemicalSystem,
    EquilibriumProblem,
    EquilibriumInverseProblem,
    isUsingOpenlibm,
)

import thermofun.PyThermoFun as thermofun

@pytest.fixture(scope="function")
def equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar():
    """
    Build a problem with 1 kg of H2O, 100 g of CO2 and 0.1 mol of NaCl
    at 60 °C and 300 bar
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhaseWithElementsOf("H2O NaCl CO2")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Halite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("CO2", 100, "g")
    problem.add("NaCl", 0.1, "mol")
    problem.setTemperature(60, "celsius")
    problem.setPressure(300, "bar")

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar():
    """
    Build a problem with H2O, H+, Na+, Cl-, HCO3-, CO2(aq), CO3-- and
    Halite at 60 °C and 300 bar
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhase(
        ["CO(aq)", "CO2(aq)", "CO3--", "Cl-", "ClO-", "ClO2-", "ClO3-", "ClO4-", "H+", "H2(aq)", "H2O(l)", "H2O2(aq)", "HCO3-", "HCl(aq)", "HClO(aq)", "HClO2(aq)", "HO2-", "Na+", "NaCl(aq)", "NaOH(aq)", "O2(aq)", "OH-"]
    ).setActivityModelDrummondCO2()
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"]).setChemicalModelSpycherPruessEnnis()
    editor.addMineralPhase("Halite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("CO2", 100, "g")
    problem.add("NaCl", 1, "mol")
    problem.setTemperature(60, "celsius")
    problem.setPressure(300, "bar")

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite():
    """
    Build a problem with H2O, Fe(OH)2, Fe(OH)3, NH3 and Magnetite
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhaseWithElementsOf("H2O Fe(OH)2 Fe(OH)3 NH3")
    editor.addGaseousPhaseWithElementsOf("NH3")
    editor.addMineralPhase("Magnetite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("Fe(OH)2", 1, "mol")
    problem.add("Fe(OH)3", 2, "mol")
    problem.add("NH3", 1, "mol")

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_problem_using_thermofun_aq17_database():
    """
    Build a problem using ThermoFun database aq17
    """

    database = thermofun.Database("databases/thermofun/aq17-thermofun.json")

    editor = ChemicalEditor(database)
    editor.setTemperatures([500.0], "celsius")
    editor.setPressures([3000.0], "bar")

    editor.addAqueousPhase([
        "Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2",
        "Ca+2", "CaCO3@", "CaCl+", "CaCl2@", "CaHCO3+", "CaHSiO3+", "CaOH+", "CaSiO3@", "K+", "KAlO2@",
        "KCl@", "KOH@", "KCO3-", "KHCO3@", "Mg+2", "MgCO3@", "MgCl+", "MgCl2@", "MgHCO3+", "MgHSiO3+",
        "MgOH+", "MgSiO3@", "Na+", "NaAl(OH)4@", "NaCO3-", "NaCl@", "NaHCO3@", "NaHSiO3@",
        "NaOH@", "HSiO3-", "SiO2@", "CO@", "CO2@", "CO3-2", "HCO3-", "CH4@", "Cl-",
        "HCl@", "H2@", "O2@", "OH-", "H+", "H2O@"])

    editor.addMineralPhase("Albite")
    editor.addMineralPhase("Andalusite")
    editor.addMineralPhase("Calcite")
    editor.addMineralPhase("Corundum")
    editor.addMineralPhase("Diopside")
    editor.addMineralPhase("Dolomite")
    editor.addMineralPhase("Enstatite")
    editor.addMineralPhase("Grossular")
    editor.addMineralPhase("Margarite")
    editor.addMineralPhase("Microcline")
    editor.addMineralPhase("Muscovite")
    editor.addMineralPhase("Pargasite-Mg")
    editor.addMineralPhase("Phlogopite")
    editor.addMineralPhase("Quartz")
    editor.addMineralPhase("Sanidine")
    editor.addMineralPhase("Sillimanite")
    editor.addMineralPhase("Zoisite")

    system = ChemicalSystem(editor)

    problem = EquilibriumProblem(system)
    problem.add("H2O",             	1000,	"g")
    problem.add("CO2",             	0.001,  "g")
    problem.add("CaCO3",           	1,	    "g")
    problem.add("MgSiO3",          	1,	    "g")
    problem.add("NaCl",            	5,      "g")
    problem.add("NaAlSi3O8",        37,	    "g")
    problem.add("KAl3Si3O10(OH)2",  13,	    "g")
    problem.add("SiO2",          	30,	    "g")
    problem.add("KAlSi3O8",        	20,	    "g")
    problem.setTemperature(500.0, "celsius")
    problem.setPressure(3000.0, "bar")

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity():
    """
    Build a problem with H, Na, Cl, Ca, Mg, C with fixed
    species amount, activity and defined pH
    """

    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C")
    editor.addGaseousPhaseWithElements("H O C")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.add("CaCl2", 2, "mmol")
    problem.add("MgCl2", 4, "mmol")
    problem.pH(3.0, "HCl")
    problem.fixSpeciesAmount("CO2(g)", 1.0, "mol")
    problem.fixSpeciesActivity("O2(g)", 0.20)

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph():
    """
    Build a problem with H, Na, Cl, Ca, Mg, C with defined pH
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhaseWithElements("H O Na Cl Ca Mg C")
    editor.addGaseousPhaseWithElements("H O C")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.add("CaCl2", 2, "mmol")
    problem.add("MgCl2", 4, "mmol")
    problem.pH(4.0, "CO2")

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts():
    """
    Build a problem with H, O, Na, Cl, Ca, C and Calcite with defined pH
    and fixed species amount
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhaseWithElements("H O Na Cl Ca C")
    editor.addMineralPhase("Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.pH(8.0, "HCl", "NaOH")
    problem.fixSpeciesAmount("Calcite", 1, "mol")

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass():
    """
    Build a problem with H2O, NaCL, CaCO3, CO2, Calcite with fixed
    species mass and amount
    """
    database = Database("supcrt98.xml")

    editor = ChemicalEditor(database)
    editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.fixSpeciesMass("Calcite", 100, "g")
    problem.fixSpeciesAmount("CO2(g)", 1.0, "mol")

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity():
    """
    Build a problem with H2O, NaCl, CaCO3, CO2 and Calcite
    with fixed values of Species Mass, Amount and alkalinity
    """
    editor = ChemicalEditor()
    editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.fixSpeciesMass("Calcite", 100, "g")
    problem.fixSpeciesAmount("CO2(g)", 1.0, "mol")
    problem.alkalinity(25.0, "meq/L", "Cl")

    return (system, problem)


@pytest.fixture(scope="function")
def equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume():
    """
    Build a problem with H2O, NaCl, CaCO3, CO2 and Calcite
    with fixed values of phase volume
    """
    editor = ChemicalEditor()
    editor.addAqueousPhaseWithElementsOf("H2O NaCl CaCO3")
    editor.addGaseousPhase(["H2O(g)", "CO2(g)"])
    editor.addMineralPhase("Calcite")

    system = ChemicalSystem(editor)

    problem = EquilibriumInverseProblem(system)
    problem.add("H2O", 1, "kg")
    problem.add("NaCl", 0.1, "mol")
    problem.fixPhaseVolume("Gaseous", 0.2, "m3", "CO2")
    problem.fixPhaseVolume("Aqueous", 0.3, "m3", "1 kg H2O; 0.1 mol NaCl")
    problem.fixPhaseVolume("Calcite", 0.5, "m3", "CaCO3")

    return (system, problem)


def _get_basename(request):
    if isUsingOpenlibm():
        import re
        return re.sub(r"[\W]", "_", request.node.name) + ".openlibm"
    return None  # use default


@pytest.fixture(scope="function")
def state_regression(num_regression, request):

    class StateRegression:
        def check(self, state, tol=None, default_tol=None, exclude=None):
            num_regression.check(
                convert_reaktoro_state_to_dict(state, exclude),
                basename=_get_basename(request),
                tolerances=tol,
                default_tolerance=default_tol,
                data_index=None,
                fill_different_shape_with_nan=True,
            )

    return StateRegression()


@pytest.fixture(scope="function")
def table_regression(num_regression, request):
    class TableRegression:
        def check(self, table, tol=None, default_tol=None):
            num_regression.check(
                convert_table_to_dict(table),
                basename=_get_basename(request),
                tolerances=tol,
                default_tolerance=default_tol,
                data_index=None,
                fill_different_shape_with_nan=True,
            )

    return TableRegression()


@pytest.fixture
def chemical_editor():
    editor = ChemicalEditor()
    editor.addAqueousPhase("H2O(l) H+ OH- HCO3- CO2(aq) CO3--".split())
    editor.addGaseousPhase("H2O(g) CO2(g)".split())
    editor.addMineralPhase("Graphite")
    return editor


@pytest.fixture
def chemical_system_adding_argon(chemical_editor):
    chemical_editor.addGaseousPhase("H2O(g) CO2(g) Ar(g)".split())
    return ChemicalSystem(chemical_editor)


@pytest.fixture
def chemical_system(chemical_editor):
    return ChemicalSystem(chemical_editor)


@pytest.fixture
def partition_with_inert_gaseous_phase(chemical_system):
    partition = Partition(chemical_system)
    partition.setInertPhases(['Gaseous'])

    return partition


@pytest.fixture
def partition_with_inert_gaseous_phase_adding_argon(chemical_system_adding_argon):
    partition = Partition(chemical_system_adding_argon)
    partition.setInertPhases(['Gaseous'])

    return partition


@pytest.fixture
def chemical_properties(chemical_system):
    # A sensible value for temperature (in K)
    T = 300

    # A sensible value for pressure (in Pa)
    P = 1e5

    # A sensible array of species amounts
    n = np.array([55, 1e-7, 1e-7, 0.1, 0.5, 0.01, 1.0, 0.001, 1.0])

    return chemical_system.properties(T, P, n)
