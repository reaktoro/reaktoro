import numpy as np
import pandas as pd
import pytest

import reaktoro

from numpy.testing import assert_allclose

is_debug_plots_on = False


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
    return df_pvtlib_result


@pytest.fixture
def composition():
    return np.array([0.5, 0.42, 0.08])


@pytest.fixture
def gaseous_species():
    return ["C1(g)", "C4(g)", "C10(g)"]


@pytest.fixture
def oil_species():
    return ["C1(liq)", "C4(liq)", "C10(liq)"]


@pytest.fixture
def custom_hydrocarbon_db(data_dir):
    return reaktoro.Database(str(data_dir / 'hydrocarbons.xml'))


@pytest.fixture
def whitson_chemical_system(gaseous_species, oil_species, custom_hydrocarbon_db):

    editor = reaktoro.ChemicalEditor(custom_hydrocarbon_db)

    eos_params = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.PengRobinson,
    )

    editor.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params)
    editor.addLiquidPhase(oil_species).setChemicalModelCubicEOS(eos_params)

    system = reaktoro.ChemicalSystem(editor)

    return system


def test_composition_results(
    composition,
    whitson_chemical_system, 
    whitson_pvtlib_results,
    num_regression
):
    temperature = 280  # degF

    system = whitson_chemical_system
    problem = reaktoro.EquilibriumProblem(system)

    problem.setTemperature(temperature, 'degF')
    problem.setElementAmount('C1', composition[0])
    problem.setElementAmount('C4', composition[1])
    problem.setElementAmount('C10', composition[2])

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

    # Retrieve pvtlib solution
    df_pvtlib_result = whitson_pvtlib_results

    df_pressures = df_pvtlib_result["P"].copy()
    df_pressures.drop_duplicates(inplace=True)

    df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.phase_id == 1]
    pvtlib_phase_fractions_gas = df_pvtlib_result_gas.F.values
    pvtlib_pressure_values_gas = df_pvtlib_result_gas.P.values

    df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.phase_id == 0]
    pvtlib_phase_fractions_liquid = df_pvtlib_result_liq.F.values
    pvtlib_pressure_values_liquid = df_pvtlib_result_liq.P.values

    pressure_values = df_pressures.values
    num_of_components = len(composition)
    num_of_phases = 2
    phase_fractions_liquid = list()
    phase_fractions_gas = list()
    composition_liq = list()
    composition_gas = list()
    pressure_values_converged = list()
    for P in pressure_values:
        problem.setPressure(P, 'psi')
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
    reaktoro_c4_gas_composition = composition_gas[:n_points_gas_pvtlib, 1]
    reaktoro_c10_gas_composition = composition_gas[:n_points_gas_pvtlib, 2]
    assert_allclose(reaktoro_c1_gas_composition, df_pvtlib_result_gas.x0, rtol=2e-3)
    assert_allclose(reaktoro_c4_gas_composition, df_pvtlib_result_gas.x1, rtol=1e-3)
    assert_allclose(reaktoro_c10_gas_composition, df_pvtlib_result_gas.x2, rtol=3e-3)

    n_points_liq_pvtlib = pvtlib_pressure_values_liquid.shape[0]
    reaktoro_c1_liq_composition = composition_liq[-n_points_liq_pvtlib:, 0]
    reaktoro_c4_liq_composition = composition_liq[-n_points_liq_pvtlib:, 1]
    reaktoro_c10_liq_composition = composition_liq[-n_points_liq_pvtlib:, 2]
    assert_allclose(reaktoro_c1_liq_composition, df_pvtlib_result_liq.x0, rtol=5e-3)
    assert_allclose(reaktoro_c4_liq_composition, df_pvtlib_result_liq.x1, rtol=2e-3)
    assert_allclose(reaktoro_c10_liq_composition, df_pvtlib_result_liq.x2, rtol=2e-3)

    molar_fractions = {
        "liquid_phase": reaktoro_phase_fractions_liquid,
        "gas_phase": reaktoro_phase_fractions_gas,
        "C1(g)": reaktoro_c1_gas_composition,
        "C4(g)": reaktoro_c4_gas_composition,
        "C10(g)": reaktoro_c10_gas_composition,
        "C1(liq)": reaktoro_c1_liq_composition,
        "C4(liq)": reaktoro_c4_liq_composition,
        "C10(liq)": reaktoro_c10_liq_composition,
    }
    num_regression.check(molar_fractions)

    if is_debug_plots_on:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 6))
        plt.plot(pressure_values_converged, phase_fractions_liquid, "-x", label="Liquid (Reaktoro)")
        plt.plot(pvtlib_pressure_values_liquid, pvtlib_phase_fractions_liquid, "-x", label="Liquid (pvtlib)")

        plt.plot(pressure_values_converged, phase_fractions_gas, "-o", label="Gas (Reaktoro)")
        plt.plot(pvtlib_pressure_values_gas, pvtlib_phase_fractions_gas, "-o", label="Gas (pvtlib)")

        plt.xlabel("Pressure [psi]")
        plt.ylabel("Phase molar fraction [mol / mol]")
        plt.title(f"Fixed T = {temperature} degF")
        plt.legend(shadow=True)

        plt.grid(True)

        # plt.savefig("reaktoro_pvtlib_using_prev.png", dpi=300)
        plt.show()

        plt.figure(figsize=(8, 6))
        plt.plot(pressure_values_converged, composition_gas[:, 0], "-x", label="C1(g) - Reaktoro")
        plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas.x0, "-x", label="C1(g) - pvtlib")
        plt.plot(pressure_values_converged, composition_gas[:, 1], "-x", label="C4(g) - Reaktoro")
        plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas.x1, "-x", label="C4(g) - pvtlib")
        plt.plot(pressure_values_converged, composition_gas[:, 2], "-x", label="C10(g) - Reaktoro")
        plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas.x2, "-x", label="C10(g) - pvtlib")

        plt.ylim([0, 1])
        plt.grid(True)

        plt.xlabel("Pressure [bar]")
        plt.ylabel("Molar fraction [mol / mol]")
        plt.title(f"Fixed T = {temperature} degC")
        plt.legend(shadow=True, ncol=3)

        # plt.savefig("reaktoro_pvtlib_compositions_whitson_gas.png", dpi=300)
        plt.show()

        plt.figure(figsize=(8, 6))
        plt.plot(pressure_values_converged, composition_liq[:, 0], "-x", label="C1(liq) - Reaktoro")
        plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq.x0, "-x", label="C1(liq) - pvtlib")
        plt.plot(pressure_values_converged, composition_liq[:, 1], "-x", label="C4(liq) - Reaktoro")
        plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq.x1, "-x", label="C4(liq) - pvtlib")
        plt.plot(pressure_values_converged, composition_liq[:, 2], "-x", label="C10(liq) - Reaktoro")
        plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq.x2, "-x", label="C10(liq) - pvtlib")

        plt.ylim([0, 1])
        plt.grid(True)

        plt.xlabel("Pressure [bar]")
        plt.ylabel("Molar fraction [mol / mol]")
        plt.title(f"Fixed T = {temperature} degC")
        plt.legend(shadow=True, ncol=3)

        # plt.savefig("reaktoro_pvtlib_compositions_whitson_liq.png", dpi=300)
        plt.show()


def test_fugacities_result(
    oil_species,
    gaseous_species,
    whitson_chemical_system,
    whitson_pvtlib_results,
    num_regression
):
    temperature = 280  # degF
    system = whitson_chemical_system

    state = reaktoro.ChemicalState(system)
    state.setTemperature(temperature, 'degF')

    # Gathering pvtlib results
    df_pvtlib_result = whitson_pvtlib_results
    df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.phase_id == 1]
    df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.phase_id == 0]

    # Properties' calculation loop for gas phase
    activities_gas = list()
    pvtlib_fugacities_gas = list()
    for _, df_row in df_pvtlib_result_gas.iterrows():
        P = df_row.P
        composition_0 = df_row.x0
        composition_1 = df_row.x1
        composition_2 = df_row.x2
        state.setPressure(P, 'psi')
        state.setSpeciesAmount(gaseous_species[0], composition_0)
        state.setSpeciesAmount(gaseous_species[1], composition_1)
        state.setSpeciesAmount(gaseous_species[2], composition_2)

        fugacity_0 = df_row.fugacity_0
        fugacity_1 = df_row.fugacity_1
        fugacity_2 = df_row.fugacity_2
        fugacities = np.array([fugacity_0, fugacity_1, fugacity_2])

        gas_properties = state.properties()
        activities = np.exp(gas_properties.lnActivities().val)
        activities_gas.append(activities)
        pvtlib_fugacities_gas.append(fugacities)

    activities_gas = np.array(activities_gas)
    pvtlib_fugacities_gas = np.array(pvtlib_fugacities_gas)

    reaktoro_c1_activity_gas = activities_gas[:, 0]
    reaktoro_c4_activity_gas = activities_gas[:, 1]
    reaktoro_c10_activity_gas = activities_gas[:, 2]
    assert_allclose(reaktoro_c1_activity_gas, pvtlib_fugacities_gas[:, 0], rtol=2e-2)
    assert_allclose(reaktoro_c4_activity_gas, pvtlib_fugacities_gas[:, 1], rtol=2e-2)
    assert_allclose(reaktoro_c10_activity_gas, pvtlib_fugacities_gas[:, 2], rtol=2e-2)

    # Properties' calculation loop for liq phase
    activities_liq = list()
    pvtlib_fugacities_liq = list()
    for _, df_row in df_pvtlib_result_liq.iterrows():
        P = df_row.P
        composition_0 = df_row.x0
        composition_1 = df_row.x1
        composition_2 = df_row.x2
        state.setPressure(P, 'psi')
        state.setSpeciesAmount(oil_species[0], composition_0)
        state.setSpeciesAmount(oil_species[1], composition_1)
        state.setSpeciesAmount(oil_species[2], composition_2)

        fugacity_0 = df_row.fugacity_0
        fugacity_1 = df_row.fugacity_1
        fugacity_2 = df_row.fugacity_2
        fugacities = np.array([fugacity_0, fugacity_1, fugacity_2])

        gas_properties = state.properties()
        activities = np.exp(gas_properties.lnActivities().val)
        activities_liq.append(activities)
        pvtlib_fugacities_liq.append(fugacities)

    activities_liq = np.array(activities_liq)
    pvtlib_fugacities_liq = np.array(pvtlib_fugacities_liq)

    reaktoro_c1_activity_liq = activities_liq[:, 3]
    reaktoro_c4_activity_liq = activities_liq[:, 4]
    reaktoro_c10_activity_liq = activities_liq[:, 5]
    assert_allclose(reaktoro_c1_activity_liq, pvtlib_fugacities_liq[:, 0], rtol=2e-2)
    assert_allclose(reaktoro_c4_activity_liq, pvtlib_fugacities_liq[:, 1], rtol=2e-2)
    assert_allclose(reaktoro_c10_activity_liq, pvtlib_fugacities_liq[:, 2], rtol=2e-2)

    components_activities = {
        "c1_activities_liq": reaktoro_c1_activity_liq,
        "c4_activities_liq": reaktoro_c4_activity_liq,
        "c10_activities_liq": reaktoro_c10_activity_liq,
        "c1_activities_gas": reaktoro_c1_activity_gas,
        "c4_activities_gas": reaktoro_c4_activity_gas,
        "c10_activities_gas": reaktoro_c10_activity_gas,
    }
    num_regression.check(components_activities)

    if is_debug_plots_on:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 6))

        pressure_values_gas = df_pvtlib_result_gas.P.values
        plt.plot(pressure_values_gas, activities_gas[:, 0], "-x", label="C1(g) - Reaktoro")
        plt.plot(pressure_values_gas, pvtlib_fugacities_gas[:, 0], "-o", label="C1(g) - pvtlib")
        plt.plot(pressure_values_gas, activities_gas[:, 1], "-x", label="C4(g) - Reaktoro")
        plt.plot(pressure_values_gas, pvtlib_fugacities_gas[:, 1], "-o", label="C4(g) - pvtlib")
        plt.plot(pressure_values_gas, activities_gas[:, 2], "-x", label="C10(g) - Reaktoro")
        plt.plot(pressure_values_gas, pvtlib_fugacities_gas[:, 2], "-o", label="C10(g) - pvtlib")

        plt.xlabel("Pressure [psi]")
        plt.ylabel("Fugacities [psi]")
        plt.title(f"Fixed T = {temperature} degF")
        plt.legend(shadow=True)

        plt.grid(True)

        # plt.savefig("reaktoro_pvtlib_fugacities_gas_whitson.png", dpi=300)
        plt.show()

        plt.figure(figsize=(8, 6))

        pressure_values_liq = df_pvtlib_result_liq.P.values
        plt.plot(pressure_values_liq, activities_liq[:, 3], "-x", label="C1(liq) - Reaktoro")
        plt.plot(pressure_values_liq, pvtlib_fugacities_liq[:, 0], "-o", label="C1(liq) - pvtlib")
        plt.plot(pressure_values_liq, activities_liq[:, 4], "-x", label="C4(liq) - Reaktoro")
        plt.plot(pressure_values_liq, pvtlib_fugacities_liq[:, 1], "-o", label="C4(liq) - pvtlib")
        plt.plot(pressure_values_liq, activities_liq[:, 5], "-x", label="C10(liq) - Reaktoro")
        plt.plot(pressure_values_liq, pvtlib_fugacities_liq[:, 2], "-o", label="C10(liq) - pvtlib")

        plt.xlabel("Pressure [psi]")
        plt.ylabel("Fugacities [psi]")
        plt.title(f"Fixed T = {temperature} degF")
        plt.legend(shadow=True)

        plt.grid(True)

        # plt.savefig("reaktoro_pvtlib_fugacities_liq_whitson.png", dpi=300)
        plt.show()
