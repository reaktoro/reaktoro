import numpy as np
import pandas as pd

import reaktoro

from numpy.testing import assert_allclose
import pytest

is_debug_plots_on = False


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
    return df_pvtlib_result


@pytest.fixture
def pedersen_pvtlib_fugacities_results(data_dir):
    df_pvtlib_result = pd.read_csv(
        data_dir / "pvtlib-pedersen63-phase-fugacities.csv"
    )
    return df_pvtlib_result


@pytest.fixture
def composition():
    return np.array([0.40, 0.60])


@pytest.fixture
def gaseous_species():
    return ["CH4(g)", "CO2(g)"]


@pytest.fixture
def oil_species():
    return ["CH4(liq)", "CO2(liq)"]


@pytest.fixture
def pedersen_chemical_system(gaseous_species, oil_species):

    temperature = -42.0  # degC

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

    editor.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params)
    editor.addLiquidPhase(oil_species).setChemicalModelCubicEOS(eos_params)

    system = reaktoro.ChemicalSystem(editor)

    return system


def test_composition_results(
    composition, 
    pedersen_chemical_system, 
    pedersen_pvtlib_composition_results,
    num_regression
):
    temperature = -42.0  # degC
    system = pedersen_chemical_system
    problem = reaktoro.EquilibriumProblem(system)

    problem.setTemperature(temperature, 'degC')
    problem.add('CH4', composition[0], "mol")
    problem.add('CO2', composition[1], "mol")
    system = problem.system()

    options = reaktoro.EquilibriumOptions()
    options.hessian = reaktoro.GibbsHessian.Exact  # required change for this case

    solver = reaktoro.EquilibriumSolver(system)
    solver.setOptions(options)

    # In this case, we should provide a better initial guess instead of the default.
    state = reaktoro.ChemicalState(system)
    state.setSpeciesAmounts(0.001)
    state.setSpeciesAmount("CH4(g)", composition[0])
    state.setSpeciesAmount("CO2(g)", composition[1])

    # Retrieve pvtlib solution
    df_pvtlib_result = pedersen_pvtlib_composition_results

    df_pressures = df_pvtlib_result["Pressure [Pa]"].copy()
    df_pressures.drop_duplicates(inplace=True)

    df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.Phase == 1]
    pvtlib_phase_fractions_gas = df_pvtlib_result_gas["Molar Fraction"].values
    pvtlib_pressure_values_gas = df_pvtlib_result_gas["Pressure [Pa]"].values
    pvtlib_pressure_values_gas = pvtlib_pressure_values_gas / 100000  # to bar

    df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.Phase == 2]
    pvtlib_phase_fractions_liquid = df_pvtlib_result_liq["Molar Fraction"].values
    pvtlib_pressure_values_liquid = df_pvtlib_result_liq["Pressure [Pa]"].values
    pvtlib_pressure_values_liquid = pvtlib_pressure_values_liquid / 100000  # to bar

    # Phase molar fractions
    temperature = problem.temperature()
    pressure_values = df_pressures.values / 100000  # to bar
    num_of_components = len(composition)
    num_of_phases = 2
    phase_fractions_liquid = list()
    phase_fractions_gas = list()
    composition_liq = list()
    composition_gas = list()
    pressure_values_converged = list()
    for P in pressure_values:
        problem.setPressure(P, 'bar')
        has_converged = solver.solve(state, problem)
        assert has_converged

        molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')

        gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
        phase_fractions_gas.append(gas_phase_molar_fraction)

        liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
        phase_fractions_liquid.append(liquid_phase_molar_fraction)

        mixture_properties = state.properties()
        x_gas = mixture_properties.moleFractions().val[0:num_of_components]
        composition_gas.append(x_gas)
        x_liq = mixture_properties.moleFractions().val[num_of_components:num_of_phases * num_of_components]
        composition_liq.append(x_liq)

        pressure_values_converged.append(P)

    phase_fractions_gas = np.array(phase_fractions_gas)
    phase_fractions_liquid = np.array(phase_fractions_liquid)
    composition_gas = np.array(composition_gas)
    composition_liq = np.array(composition_liq)
    pressure_values_converged = np.array(pressure_values_converged)

    zero_threshold = 1e-5

    # Phase fractions
    reaktoro_phase_fractions_liquid = phase_fractions_liquid[phase_fractions_liquid > zero_threshold]
    reaktoro_phase_fractions_gas = phase_fractions_gas[phase_fractions_gas > zero_threshold]

    assert_allclose(reaktoro_phase_fractions_liquid, pvtlib_phase_fractions_liquid, atol=1e-2)
    assert_allclose(reaktoro_phase_fractions_gas, pvtlib_phase_fractions_gas, atol=1e-2)

    # Component molar fractions
    n_points_gas_pvtlib = pvtlib_pressure_values_gas.shape[0]
    reaktoro_c1_gas_composition = composition_gas[:n_points_gas_pvtlib, 0]
    reaktoro_co2_gas_composition = composition_gas[:n_points_gas_pvtlib, 1]
    assert_allclose(reaktoro_c1_gas_composition, df_pvtlib_result_gas["C1"], rtol=2e-3)
    assert_allclose(reaktoro_co2_gas_composition, df_pvtlib_result_gas["CO2"], rtol=1e-3)

    n_points_liq_pvtlib = pvtlib_pressure_values_liquid.shape[0]
    reaktoro_c1_liq_composition = composition_liq[-n_points_liq_pvtlib:, 0]
    reaktoro_co2_liq_composition = composition_liq[-n_points_liq_pvtlib:, 1]
    assert_allclose(reaktoro_c1_liq_composition, df_pvtlib_result_liq["C1"], rtol=5e-3)
    assert_allclose(reaktoro_co2_liq_composition, df_pvtlib_result_liq["CO2"], rtol=2e-3)

    molar_fractions = {
        "liquid_phase": reaktoro_phase_fractions_liquid,
        "gas_phase": reaktoro_phase_fractions_gas,
        "C1(g)": reaktoro_c1_gas_composition,
        "CO2(g)": reaktoro_co2_gas_composition,
        "C1(liq)": reaktoro_c1_liq_composition,
        "CO2(liq)": reaktoro_co2_liq_composition,
    }
    num_regression.check(molar_fractions)

    if is_debug_plots_on:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8, 6))
        plt.plot(pressure_values_converged, phase_fractions_liquid, "-x", label="Liquid (Reaktoro)")
        plt.plot(pvtlib_pressure_values_liquid, pvtlib_phase_fractions_liquid, "-x", label="Liquid (pvtlib)")

        plt.plot(pressure_values_converged, phase_fractions_gas, "-o", label="Gas (Reaktoro)")
        plt.plot(pvtlib_pressure_values_gas, pvtlib_phase_fractions_gas, "-o", label="Gas (pvtlib)")

        plt.xlabel("Pressure [bar]")
        plt.ylabel("Phase molar fraction [mol / mol]")
        plt.title(f"Fixed T = {temperature} degC")
        plt.legend(shadow=True)

        plt.grid(True)

        # plt.savefig("reaktoro_pvtlib_using_prev_pedersen.png", dpi=300)
        plt.show()

        plt.figure(figsize=(8, 6))
        plt.plot(pvtlib_pressure_values_gas, reaktoro_c1_gas_composition, "-x", label="C1(g) - Reaktoro")
        plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas["C1"], "-x", label="C1(g) - pvtlib")
        plt.plot(pvtlib_pressure_values_gas, reaktoro_co2_gas_composition, "-x", label="CO2(g) - Reaktoro")
        plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas["CO2"], "-x", label="CO2(g) - pvtlib")

        plt.ylim([0, 1])
        plt.grid(True)

        plt.xlabel("Pressure [bar]")
        plt.ylabel("Molar fraction [mol / mol]")
        plt.title(f"Fixed T = {temperature} degC")
        plt.legend(shadow=True, ncol=2)

        # plt.savefig("reaktoro_pvtlib_compositions_pedersen_gas.png", dpi=300)
        plt.show()

        plt.figure(figsize=(8, 6))
        plt.plot(pvtlib_pressure_values_liquid, reaktoro_c1_liq_composition, "-x", label="C1(liq) - Reaktoro")
        plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq["C1"], "-x", label="C1(liq) - pvtlib")
        plt.plot(pvtlib_pressure_values_liquid, reaktoro_co2_liq_composition, "-x", label="CO2(liq) - Reaktoro")
        plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq["CO2"], "-x", label="CO2(liq) - pvtlib")

        plt.ylim([0, 1])
        plt.grid(True)

        plt.xlabel("Pressure [bar]")
        plt.ylabel("Molar fraction [mol / mol]")
        plt.title(f"Fixed T = {temperature} degC")
        plt.legend(shadow=True, ncol=2)

        # plt.savefig("reaktoro_pvtlib_compositions_pedersen_liq.png", dpi=300)
        plt.show()


def test_fugacities_results(
    oil_species,
    gaseous_species,
    pedersen_chemical_system, 
    pedersen_pvtlib_composition_results,
    pedersen_pvtlib_fugacities_results,
    num_regression
):
    temperature = -42.0  # degC
    system = pedersen_chemical_system
    state = reaktoro.ChemicalState(system)
    state.setTemperature(temperature, 'degC')

    # Gathering pvtlib results
    df_pvtlib_fugacities = pedersen_pvtlib_fugacities_results
    df_pvtlib_compositions = pedersen_pvtlib_composition_results

    # Normalizing fugacities
    P_1_atm_as_Pa = 101325.0
    df_pvtlib_fugacities["C1"] = df_pvtlib_fugacities["C1"] / P_1_atm_as_Pa
    df_pvtlib_fugacities["CO2"] = df_pvtlib_fugacities["CO2"] / P_1_atm_as_Pa

    # Retrieve results by phase
    df_pvtlib_fugacities_gas = df_pvtlib_fugacities[df_pvtlib_fugacities.Phase == 1]
    df_pvtlib_fugacities_gas.reset_index(drop=True, inplace=True)
    df_pvtlib_fugacities_liq = df_pvtlib_fugacities[df_pvtlib_fugacities.Phase == 2]
    df_pvtlib_fugacities_liq.reset_index(drop=True, inplace=True)
    df_pvtlib_composition_gas = df_pvtlib_compositions[df_pvtlib_compositions.Phase == 1]
    df_pvtlib_composition_gas.reset_index(drop=True, inplace=True)
    df_pvtlib_composition_liq = df_pvtlib_compositions[df_pvtlib_compositions.Phase == 2]
    df_pvtlib_composition_liq.reset_index(drop=True, inplace=True)

    # Properties' calculation loop for gas phase
    activities_gas = list()
    pvtlib_fugacities_gas = list()
    for index, df_row in df_pvtlib_composition_gas.iterrows():
        P = df_row["Pressure [Pa]"]
        composition_0 = df_row["C1"]
        composition_1 = df_row["CO2"]
        state.setPressure(P, 'Pa')
        state.setSpeciesAmount(gaseous_species[0], composition_0)
        state.setSpeciesAmount(gaseous_species[1], composition_1)

        df_row_fugacities = df_pvtlib_fugacities_gas.iloc[index, :]
        fugacity_0 = df_row_fugacities["C1"]
        fugacity_1 = df_row_fugacities["CO2"]
        fugacities = np.array([fugacity_0, fugacity_1])

        gas_properties = state.properties()
        activities = np.exp(gas_properties.lnActivities().val)
        activities_gas.append(activities)
        pvtlib_fugacities_gas.append(fugacities)

    activities_gas = np.array(activities_gas)
    pvtlib_fugacities_gas = np.array(pvtlib_fugacities_gas)

    reaktoro_c1_activity_gas = activities_gas[:, 0]
    reaktoro_co2_activity_gas = activities_gas[:, 1]
    assert_allclose(reaktoro_c1_activity_gas, pvtlib_fugacities_gas[:, 0], rtol=2e-2)
    assert_allclose(reaktoro_co2_activity_gas, pvtlib_fugacities_gas[:, 1], rtol=2e-2)

    # Properties' calculation loop for gas phase
    activities_liq = list()
    pvtlib_fugacities_liq = list()
    for index, df_row in df_pvtlib_composition_liq.iterrows():
        P = df_row["Pressure [Pa]"]
        composition_0 = df_row["C1"]
        composition_1 = df_row["CO2"]
        state.setPressure(P, 'Pa')
        state.setSpeciesAmount(oil_species[0], composition_0)
        state.setSpeciesAmount(oil_species[1], composition_1)

        df_row_fugacities = df_pvtlib_fugacities_liq.iloc[index, :]
        fugacity_0 = df_row_fugacities["C1"]
        fugacity_1 = df_row_fugacities["CO2"]
        fugacities = np.array([fugacity_0, fugacity_1])

        liq_properties = state.properties()
        activities = np.exp(liq_properties.lnActivities().val)
        activities_liq.append(activities)
        pvtlib_fugacities_liq.append(fugacities)

    activities_liq = np.array(activities_liq)
    pvtlib_fugacities_liq = np.array(pvtlib_fugacities_liq)

    reaktoro_c1_activity_liq = activities_liq[:, 2]
    reaktoro_co2_activity_liq = activities_liq[:, 3]
    assert_allclose(reaktoro_c1_activity_liq, pvtlib_fugacities_liq[:, 0], rtol=2e-2)
    assert_allclose(reaktoro_co2_activity_liq, pvtlib_fugacities_liq[:, 1], rtol=2e-2)

    components_activities = {
        "c1_activities_liq": reaktoro_c1_activity_liq,
        "co2_activities_liq": reaktoro_co2_activity_liq,
        "c1_activities_gas": reaktoro_c1_activity_gas,
        "co2_activities_gas": reaktoro_co2_activity_gas,
    }
    num_regression.check(components_activities)

    if is_debug_plots_on:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 6))

        pressure_values_gas = df_pvtlib_composition_gas["Pressure [Pa]"].values / 100000  # to bar
        plt.plot(pressure_values_gas, activities_gas[:, 0], "-x", label="C1(g) - Reaktoro")
        plt.plot(pressure_values_gas, pvtlib_fugacities_gas[:, 0], "-o", label="C1(g) - pvtlib")
        plt.plot(pressure_values_gas, activities_gas[:, 1], "-x", label="CO2(g) - Reaktoro")
        plt.plot(pressure_values_gas, pvtlib_fugacities_gas[:, 1], "-o", label="CO2(g) - pvtlib")

        plt.xlabel("Pressure [bar]")
        plt.ylabel("Fugacities [bar]")
        plt.title(f"Fixed T = {temperature} degC")
        plt.legend(shadow=True)

        plt.grid(True)

        # plt.savefig("reaktoro_pvtlib_fugacities_gas_pedersen.png", dpi=300)
        plt.show()

        plt.figure(figsize=(8, 6))

        pressure_values_liq = df_pvtlib_composition_liq["Pressure [Pa]"].values / 100000  # to bar
        plt.plot(pressure_values_liq, activities_liq[:, 2], "-x", label="C1(liq) - Reaktoro")
        plt.plot(pressure_values_liq, pvtlib_fugacities_liq[:, 0], "-o", label="C1(liq) - pvtlib")
        plt.plot(pressure_values_liq, activities_liq[:, 3], "-x", label="CO2(liq) - Reaktoro")
        plt.plot(pressure_values_liq, pvtlib_fugacities_liq[:, 1], "-o", label="CO2(liq) - pvtlib")

        plt.xlabel("Pressure [bar]")
        plt.ylabel("Fugacities [bar]")
        plt.title(f"Fixed T = {temperature} degC")
        plt.legend(shadow=True)

        plt.grid(True)

        # plt.savefig("reaktoro_pvtlib_fugacities_liq_pedersen.png", dpi=300)
        plt.show()
