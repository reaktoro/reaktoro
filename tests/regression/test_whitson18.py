import numpy as np
import pandas as pd
import pytest

import reaktoro

from numpy.testing import assert_allclose


@pytest.fixture
def data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parents[1] / "data"


@pytest.fixture
def whitson_pvtlib_results(data_dir):
    df_pvtlib_result = pd.read_csv(
        data_dir / "pvtlib_table_whitson18_fixed_T.csv",
    )

    whitson_results = dict()

    df_pressures = df_pvtlib_result["P"].copy()
    df_pressures.drop_duplicates(inplace=True)

    whitson_results["Simulated Pressure [psi]"] = np.array(df_pressures.values)

    df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.phase_id == 1]
    whitson_results["Phase Fraction - gas"] = np.array(df_pvtlib_result_gas.F.values)
    whitson_results["Pressure - gas [psi]"] = np.array(df_pvtlib_result_gas.P.values)
    whitson_results["C1(g)"] = np.array(df_pvtlib_result_gas.x0)
    whitson_results["C4(g)"] = np.array(df_pvtlib_result_gas.x1)
    whitson_results["C10(g)"] = np.array(df_pvtlib_result_gas.x2)
    whitson_results["fugacities C1(g)"] = np.array(df_pvtlib_result_gas.fugacity_0)
    whitson_results["fugacities C4(g)"] = np.array(df_pvtlib_result_gas.fugacity_1)
    whitson_results["fugacities C10(g)"] = np.array(df_pvtlib_result_gas.fugacity_2)

    df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.phase_id == 0]
    whitson_results["Phase Fraction - liq"] = np.array(df_pvtlib_result_liq.F.values)
    whitson_results["Pressure - liq [psi]"] = np.array(df_pvtlib_result_liq.P.values)
    whitson_results["C1(liq)"] = np.array(df_pvtlib_result_liq.x0)
    whitson_results["C4(liq)"] = np.array(df_pvtlib_result_liq.x1)
    whitson_results["C10(liq)"] = np.array(df_pvtlib_result_liq.x2)
    whitson_results["fugacities C1(liq)"] = np.array(df_pvtlib_result_liq.fugacity_0)
    whitson_results["fugacities C4(liq)"] = np.array(df_pvtlib_result_liq.fugacity_1)
    whitson_results["fugacities C10(liq)"] = np.array(df_pvtlib_result_liq.fugacity_2)

    return whitson_results


@pytest.fixture
def composition():
    return np.array([0.5, 0.42, 0.08])


@pytest.fixture
def custom_hydrocarbon_db(data_dir):
    return reaktoro.Database(str(data_dir / 'hydrocarbons.xml'))


@pytest.fixture
def whitson_chemical_system(custom_hydrocarbon_db):

    editor = reaktoro.ChemicalEditor(custom_hydrocarbon_db)

    eos_params = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.PengRobinson,
    )

    editor.addGaseousPhase(["C1(g)", "C4(g)", "C10(g)"]).setChemicalModelCubicEOS(eos_params)
    editor.addLiquidPhase(["C1(liq)", "C4(liq)", "C10(liq)"]).setChemicalModelCubicEOS(eos_params)

    system = reaktoro.ChemicalSystem(editor)

    return system


def test_composition_results(
    composition,
    whitson_chemical_system, 
    whitson_pvtlib_results,
    num_regression
):
    system = whitson_chemical_system

    options = reaktoro.EquilibriumOptions()
    options.hessian = reaktoro.GibbsHessian.Exact  # required change for this case

    solver = reaktoro.EquilibriumSolver(system)
    solver.setOptions(options)

    # In this case, we should provide a better initial guess instead of the default.
    state = reaktoro.ChemicalState(system)
    state.setSpeciesAmounts(0.001)  # start will all having 0.001 moles
    state.setSpeciesAmount("C1(g)", composition[0])  # overwrite amount of C1(g) (same below)
    state.setSpeciesAmount("C4(g)", composition[1])  
    state.setSpeciesAmount("C10(g)", composition[2])

    temperature = 280  # degF

    problem = reaktoro.EquilibriumProblem(system)
    problem.setTemperature(temperature, 'degF')
    problem.setElementAmount('C1', composition[0])
    problem.setElementAmount('C4', composition[1])
    problem.setElementAmount('C10', composition[2])

    phase_fractions_liquid = list()
    phase_fractions_gas = list()
    reaktoro_composition_liq = {
        'C1(liq)': list(),
        'C4(liq)': list(),
        'C10(liq)': list(),
    }
    reaktoro_composition_gas = {
        'C1(g)': list(),
        'C4(g)': list(),
        'C10(g)': list(),
    }
    pressure_values = whitson_pvtlib_results["Simulated Pressure [psi]"]
    for pressure in pressure_values:
        problem.setPressure(pressure, 'psi')
        has_converged = solver.solve(state, problem)
        assert has_converged

        molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')

        gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
        phase_fractions_gas.append(gas_phase_molar_fraction)

        liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
        phase_fractions_liquid.append(liquid_phase_molar_fraction)

        mixture_properties = state.properties()
        if (pressure in whitson_pvtlib_results["Pressure - gas [psi]"]):
            reaktoro_composition_gas["C1(g)"].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("C1(g)")])
            reaktoro_composition_gas["C4(g)"].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("C4(g)")])
            reaktoro_composition_gas["C10(g)"].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("C10(g)")])

        if (pressure in whitson_pvtlib_results["Pressure - liq [psi]"]):
            reaktoro_composition_liq["C1(liq)"].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("C1(liq)")])
            reaktoro_composition_liq["C4(liq)"].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("C4(liq)")])
            reaktoro_composition_liq["C10(liq)"].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("C10(liq)")])

    phase_fractions_gas = np.array(phase_fractions_gas)
    phase_fractions_liquid = np.array(phase_fractions_liquid)

    zero_threshold = 1e-5

    # Component Phase fractions
    reaktoro_phase_fractions_liquid = phase_fractions_liquid[phase_fractions_liquid > zero_threshold]
    reaktoro_phase_fractions_gas = phase_fractions_gas[phase_fractions_gas > zero_threshold]

    assert_allclose(reaktoro_phase_fractions_liquid, whitson_pvtlib_results["Phase Fraction - liq"], atol=1e-2)
    assert_allclose(reaktoro_phase_fractions_gas, whitson_pvtlib_results["Phase Fraction - gas"], atol=1e-2)

    # Compare phase fractions with results computed by pvtlib
    assert_allclose(reaktoro_composition_gas["C1(g)"], whitson_pvtlib_results["C1(g)"], rtol=2e-3)
    assert_allclose(reaktoro_composition_gas["C4(g)"], whitson_pvtlib_results["C4(g)"], rtol=1e-3)
    assert_allclose(reaktoro_composition_gas["C10(g)"], whitson_pvtlib_results["C10(g)"], rtol=3e-3)

    assert_allclose(reaktoro_composition_liq["C1(liq)"], whitson_pvtlib_results["C1(liq)"], rtol=5e-3)
    assert_allclose(reaktoro_composition_liq["C4(liq)"], whitson_pvtlib_results["C4(liq)"], rtol=2e-3)
    assert_allclose(reaktoro_composition_liq["C10(liq)"], whitson_pvtlib_results["C10(liq)"], rtol=2e-3)

    molar_fractions = {
        "liquid_phase": reaktoro_phase_fractions_liquid,
        "gas_phase": reaktoro_phase_fractions_gas,
        "C1(g)": reaktoro_composition_gas["C1(g)"],
        "C4(g)": reaktoro_composition_gas["C4(g)"],
        "C10(g)": reaktoro_composition_gas["C10(g)"],
        "C1(liq)": reaktoro_composition_liq["C1(liq)"],
        "C4(liq)": reaktoro_composition_liq["C4(liq)"],
        "C10(liq)": reaktoro_composition_liq["C10(liq)"],
    }
    num_regression.check(molar_fractions)

def test_fugacities_result(
    whitson_chemical_system,
    whitson_pvtlib_results,
    num_regression
):
    temperature = 280  # degF
    system = whitson_chemical_system

    state = reaktoro.ChemicalState(system)
    state.setTemperature(temperature, 'degF')

    # Properties' calculation loop for gas phase
    reaktoro_activities_gas = {
        'C1(g)': list(),
        'C4(g)': list(),
        'C10(g)': list(),
    }
    for index in range(len(whitson_pvtlib_results["Pressure - gas [psi]"])):
        pressure = whitson_pvtlib_results["Pressure - gas [psi]"][index]
        c1_g = whitson_pvtlib_results["C1(g)"][index]
        c4_g = whitson_pvtlib_results["C4(g)"][index]
        c10_g = whitson_pvtlib_results["C10(g)"][index]
        state.setPressure(pressure, 'psi')
        state.setSpeciesAmount("C1(g)", c1_g)
        state.setSpeciesAmount("C4(g)", c4_g)
        state.setSpeciesAmount("C10(g)", c10_g)

        gas_properties = state.properties()
        activities = np.exp(gas_properties.lnActivities().val)
        reaktoro_activities_gas['C1(g)'].append(activities[0])
        reaktoro_activities_gas['C4(g)'].append(activities[1])
        reaktoro_activities_gas['C10(g)'].append(activities[2])

    assert_allclose(reaktoro_activities_gas['C1(g)'], whitson_pvtlib_results["fugacities C1(g)"], rtol=2e-2)
    assert_allclose(reaktoro_activities_gas['C4(g)'], whitson_pvtlib_results["fugacities C4(g)"], rtol=2e-2)
    assert_allclose(reaktoro_activities_gas['C10(g)'], whitson_pvtlib_results["fugacities C10(g)"], rtol=2e-2)

    # Properties' calculation loop for liq phase
    reaktoro_activities_liq = {
        'C1(liq)': list(),
        'C4(liq)': list(),
        'C10(liq)': list(),
    }
    for index in range(len(whitson_pvtlib_results["Pressure - liq [psi]"])):
        pressure = whitson_pvtlib_results["Pressure - liq [psi]"][index]
        c1_liq = whitson_pvtlib_results["C1(liq)"][index]
        c4_liq = whitson_pvtlib_results["C4(liq)"][index]
        c10_liq = whitson_pvtlib_results["C10(liq)"][index]
        state.setPressure(pressure, 'psi')
        state.setSpeciesAmount('C1(liq)', c1_liq)
        state.setSpeciesAmount('C4(liq)', c4_liq)
        state.setSpeciesAmount('C10(liq)', c10_liq)

        liq_properties = state.properties()
        activities = np.exp(liq_properties.lnActivities().val)
        reaktoro_activities_liq['C1(liq)'].append(activities[3])
        reaktoro_activities_liq['C4(liq)'].append(activities[4])
        reaktoro_activities_liq['C10(liq)'].append(activities[5])

    assert_allclose(reaktoro_activities_liq['C1(liq)'], whitson_pvtlib_results["fugacities C1(liq)"], rtol=2e-2)
    assert_allclose(reaktoro_activities_liq['C4(liq)'], whitson_pvtlib_results["fugacities C4(liq)"], rtol=2e-2)
    assert_allclose(reaktoro_activities_liq['C10(liq)'], whitson_pvtlib_results["fugacities C10(liq)"], rtol=2e-2)

    components_activities = {
        "c1_activities_liq": reaktoro_activities_liq['C1(liq)'],
        "c4_activities_liq": reaktoro_activities_liq['C4(liq)'],
        "c10_activities_liq": reaktoro_activities_liq['C10(liq)'],
        "c1_activities_gas": reaktoro_activities_gas['C1(g)'],
        "c4_activities_gas": reaktoro_activities_gas['C4(g)'],
        "c10_activities_gas": reaktoro_activities_gas['C10(g)'],
    }
    num_regression.check(components_activities)
