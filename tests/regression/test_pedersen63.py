import numpy as np
import pandas as pd

import reaktoro

from numpy.testing import assert_allclose
import pytest


@pytest.fixture
def data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parents[1] / "data"


@pytest.fixture
def pedersen_pvtlib_composition_results(data_dir):
    df_pvtlib_result = pd.read_csv(
        data_dir / "pvtlib-pedersen63-phase-compositions.csv",
    )
    pedersen_results = dict()

    df_pressures = df_pvtlib_result["Pressure [Pa]"].copy()
    df_pressures.drop_duplicates(inplace=True)

    pedersen_results["Simulated Pressure [bar]"] = np.array(df_pressures.values/ 100000)

    df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.Phase == 1]
    pedersen_results["Molar Fraction - gas"] = np.array(df_pvtlib_result_gas["Molar Fraction"].values)
    pedersen_results["Pressure - gas [bar]"] = np.array(df_pvtlib_result_gas["Pressure [Pa]"].values / 100000)
    pedersen_results["CH4(g)"] = np.array(df_pvtlib_result_gas["C1"].values)
    pedersen_results["CO2(g)"] = np.array(df_pvtlib_result_gas["CO2"].values)

    df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.Phase == 2]
    pedersen_results["Molar Fraction - liq"] = np.array(df_pvtlib_result_liq["Molar Fraction"].values)
    pedersen_results["Pressure - liq [bar]"] = np.array(df_pvtlib_result_liq["Pressure [Pa]"].values / 100000)
    pedersen_results["CH4(liq)"] = np.array(df_pvtlib_result_liq["C1"].values)
    pedersen_results["CO2(liq)"] = np.array(df_pvtlib_result_liq["CO2"].values)

    return pedersen_results


@pytest.fixture
def pedersen_pvtlib_fugacities_results(data_dir):
    df_pvtlib_result = pd.read_csv(
        data_dir / "pvtlib-pedersen63-phase-fugacities.csv"
    )
    pedersen_results = dict()

    df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.Phase == 1]

    pedersen_results["fugacities CH4(g)"] = np.array(df_pvtlib_result_gas["C1"].values)
    pedersen_results["fugacities CO2(g)"] = np.array(df_pvtlib_result_gas["CO2"].values)

    df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.Phase == 2]

    pedersen_results["fugacities CH4(liq)"] = np.array(df_pvtlib_result_liq["C1"].values)
    pedersen_results["fugacities CO2(liq)"] = np.array(df_pvtlib_result_liq["CO2"].values)

    return pedersen_results


@pytest.fixture
def composition():
    return np.array([0.40, 0.60])


@pytest.fixture
def pedersen_chemical_system():
    db = reaktoro.Database('supcrt98.xml')
    editor = reaktoro.ChemicalEditor(db)

    def calculate_bips(T):
        k = np.zeros((2, 2))
        k[0, 1] = 0.12
        k[1, 0] = 0.12
        bips = reaktoro.BinaryInteractionParams(k)
        return bips

    eos_params = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.SoaveRedlichKwong,
        binary_interaction_values=calculate_bips
    )

    editor.addGaseousPhase(["CH4(g)", "CO2(g)"]).setChemicalModelCubicEOS(eos_params)
    editor.addLiquidPhase(["CH4(liq)", "CO2(liq)"]).setChemicalModelCubicEOS(eos_params)

    system = reaktoro.ChemicalSystem(editor)

    return system


def test_composition_results(
    composition, 
    pedersen_chemical_system, 
    pedersen_pvtlib_composition_results,
    num_regression
):
    system = pedersen_chemical_system

    options = reaktoro.EquilibriumOptions()
    options.hessian = reaktoro.GibbsHessian.Exact  # required change for this case

    solver = reaktoro.EquilibriumSolver(system)
    solver.setOptions(options)

    # In this case, we should provide a better initial guess instead of the default.
    state = reaktoro.ChemicalState(system)
    state.setSpeciesAmounts(0.001)
    state.setSpeciesAmount("CH4(g)", composition[0])
    state.setSpeciesAmount("CO2(g)", composition[1])

    temperature = -42.0  # degC
    problem = reaktoro.EquilibriumProblem(system)
    problem.addState(state)
    problem.setTemperature(temperature, 'degC')

    # Phase molar fractions
    phase_fractions_liquid = list()
    phase_fractions_gas = list()
    reaktoro_composition_liq = {
        'CH4(liq)': list(),
        'CO2(liq)': list(),
    }
    reaktoro_composition_gas = {
        'CH4(g)': list(),
        'CO2(g)': list(),
    }
    pressure_values = pedersen_pvtlib_composition_results["Simulated Pressure [bar]"]
    for pressure in pressure_values:
        problem.setPressure(pressure, 'bar')
        has_converged = solver.solve(state, problem)
        assert has_converged

        molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')

        gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
        phase_fractions_gas.append(gas_phase_molar_fraction)

        liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
        phase_fractions_liquid.append(liquid_phase_molar_fraction)

        mixture_properties = state.properties()
        if (pressure in pedersen_pvtlib_composition_results["Pressure - liq [bar]"]):
            reaktoro_composition_liq['CH4(liq)'].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("CH4(liq)")])
            reaktoro_composition_liq['CO2(liq)'].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("CO2(liq)")])
        if (pressure in pedersen_pvtlib_composition_results["Pressure - gas [bar]"]):
            reaktoro_composition_gas['CH4(g)'].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("CH4(g)")])
            reaktoro_composition_gas['CO2(g)'].append(
                mixture_properties.moleFractions().val[state.system().indexSpecies("CO2(g)")])

    phase_fractions_gas = np.array(phase_fractions_gas)
    phase_fractions_liquid = np.array(phase_fractions_liquid)

    zero_threshold = 1e-5

    # Component Phase fractions
    reaktoro_phase_fractions_liquid = phase_fractions_liquid[phase_fractions_liquid > zero_threshold]
    reaktoro_phase_fractions_gas = phase_fractions_gas[phase_fractions_gas > zero_threshold]

    # Compare phase fractions with results computed by pvtlib
    assert_allclose(reaktoro_phase_fractions_liquid, pedersen_pvtlib_composition_results["Molar Fraction - liq"], atol=1e-2)
    assert_allclose(reaktoro_phase_fractions_gas, pedersen_pvtlib_composition_results["Molar Fraction - gas"], atol=1e-2)

    # Compare component molar fractions with results computed by pvtlib
    assert_allclose(reaktoro_composition_gas['CH4(g)'], pedersen_pvtlib_composition_results['CH4(g)'], rtol=2e-3)
    assert_allclose(reaktoro_composition_gas['CO2(g)'], pedersen_pvtlib_composition_results['CO2(g)'], rtol=1e-3)
    assert_allclose(reaktoro_composition_liq['CH4(liq)'], pedersen_pvtlib_composition_results['CH4(liq)'], rtol=5e-3)
    assert_allclose(reaktoro_composition_liq['CO2(liq)'], pedersen_pvtlib_composition_results['CO2(liq)'], rtol=2e-3)

    molar_fractions = {
        "liquid_phase": reaktoro_phase_fractions_liquid,
        "gas_phase": reaktoro_phase_fractions_gas,
        "C1(g)": reaktoro_composition_gas['CH4(g)'],
        "CO2(g)": reaktoro_composition_gas['CO2(g)'],
        "C1(liq)": reaktoro_composition_liq['CH4(liq)'],
        "CO2(liq)": reaktoro_composition_liq['CO2(liq)'],
    }
    num_regression.check(molar_fractions)


def test_fugacities_results(
    pedersen_chemical_system,
    pedersen_pvtlib_composition_results,
    pedersen_pvtlib_fugacities_results,
    num_regression
):
    temperature = -42.0  # degC
    system = pedersen_chemical_system
    state = reaktoro.ChemicalState(system)
    state.setTemperature(temperature, 'degC')

    # Properties' calculation loop for gas phase
    activities_ch4_g = list()
    activities_co2_g = list()
    for index in range(len(pedersen_pvtlib_composition_results["Pressure - gas [bar]"])):
        pressure = pedersen_pvtlib_composition_results["Pressure - gas [bar]"][index]
        ch4_g = pedersen_pvtlib_composition_results["CH4(g)"][index]
        co2_g = pedersen_pvtlib_composition_results["CO2(g)"][index]
        state.setPressure(pressure, 'bar')
        state.setSpeciesAmount("CH4(g)", ch4_g)
        state.setSpeciesAmount("CO2(g)", co2_g)

        gas_properties = state.properties()
        activities = np.exp(gas_properties.lnActivities().val)
        activities_ch4_g.append(activities[0])
        activities_co2_g.append(activities[1])

    P_1_atm_as_Pa = 101325.0
    assert_allclose(activities_ch4_g, pedersen_pvtlib_fugacities_results["fugacities CH4(g)"]/P_1_atm_as_Pa, rtol=2e-2)
    assert_allclose(activities_co2_g, pedersen_pvtlib_fugacities_results["fugacities CO2(g)"]/P_1_atm_as_Pa, rtol=2e-2)

    # Properties' calculation loop for gas phase
    activities_ch4_liq = list()
    activities_co2_liq = list()
    for index in range(len(pedersen_pvtlib_composition_results["Pressure - liq [bar]"])):
        pressure = pedersen_pvtlib_composition_results["Pressure - liq [bar]"][index]
        ch4_liq = pedersen_pvtlib_composition_results["CH4(liq)"][index]
        co2_liq = pedersen_pvtlib_composition_results["CO2(liq)"][index]
        state.setPressure(pressure, 'bar')
        state.setSpeciesAmount("CH4(liq)", ch4_liq)
        state.setSpeciesAmount("CO2(liq)", co2_liq)

        liq_properties = state.properties()
        activities = np.exp(liq_properties.lnActivities().val)
        activities_ch4_liq.append(activities[2])
        activities_co2_liq.append(activities[3])

    assert_allclose(activities_ch4_liq, pedersen_pvtlib_fugacities_results["fugacities CH4(liq)"]/P_1_atm_as_Pa, rtol=2e-2)
    assert_allclose(activities_co2_liq, pedersen_pvtlib_fugacities_results["fugacities CO2(liq)"]/P_1_atm_as_Pa, rtol=2e-2)

    components_activities = {
        "ch4_liq_activities": activities_ch4_liq,
        "co2_liq_activities": activities_co2_liq,
        "ch4_gas_activities": activities_ch4_g,
        "co2_gas_activities": activities_co2_g,
    }
    num_regression.check(components_activities)
