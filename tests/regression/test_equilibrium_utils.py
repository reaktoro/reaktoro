import numpy as np
import pytest
from reaktoro import ChemicalState, EquilibriumOptions, Partition, equilibrate


# DO NOT TRY TO PUT ALL EQUILIBRIUM TEST IN A SINGUE TEST
# if one one of then fail the test will stop and you won't test all
@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume"
            )
        ),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite at 60 °C and 300 bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved at 60 °C and 300 bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrate_overload_1(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of a problem using
    equilibrate(const EquilibriumProblem& problem)
    and
    equilibrate(const EquilibriumInverseProblem& problem)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """

    (system, problem) = setup

    equilibriumState = equilibrate(problem)

    state_regression.check(equilibriumState, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume"
            )
        ),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite at 60 °C and 300 bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved at 60 °C and 300 bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrate_overload_2(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of a problem using
    equilibrate(const EquilibriumProblem& problem, const EquilibriumOptions& options)
    and
    equilibrate(const EquilibriumInverseProblem& problem, const EquilibriumOptions& options)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    options = EquilibriumOptions()

    equilibriumState = equilibrate(problem, options)

    state_regression.check(equilibriumState, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume"
            )
        ),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite at 60 °C and 300 bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved at 60 °C and 300 bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrate_overload_3(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of a problem using
    equilibrate(ChemicalState& state)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """

    (system, problem) = setup

    # find a state based on the equilibrium setup
    equilibriumState = equilibrate(problem)

    # change pressure and temperature
    equilibriumState.setTemperature(problem.temperature() + 5)
    equilibriumState.setPressure(problem.pressure() + 5)

    # compute equilibrium for state with new temperature and pressure
    equilibriumResult = equilibrate(equilibriumState)

    state_regression.check(equilibriumState, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume"
            )
        ),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite at 60 °C and 300 bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved at 60 °C and 300 bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrate_overload_4(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of a problem using
    equilibrate(ChemicalState& state, const Partition& partition)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    # find a state based on the equilibrium setup
    state = equilibrate(problem)

    # change pressure and temperature
    state.setTemperature(problem.temperature() + 5)
    state.setPressure(problem.pressure() + 5)

    partition = Partition(system)

    # compute equilibrium for state with new temperature and pressure
    equilibriumResult = equilibrate(state, partition)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume"
            )
        ),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite at 60 °C and 300 bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved at 60 °C and 300 bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrate_overload_5(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of a problem using
    equilibrate(ChemicalState& state, const EquilibriumOptions& options)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    # Find a state based on the equilibrium setup
    state = equilibrate(problem)

    # Change pressure and temperature
    state.setTemperature(problem.temperature() + 5)
    state.setPressure(problem.pressure() + 5)

    options = EquilibriumOptions()

    # Compute equilibrium for state with new temperature and pressure
    equilibriumResult = equilibrate(state, options)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume"
            )
        ),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite at 60 °C and 300 bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved at 60 °C and 300 bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrate_overload_6(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of a problem using
    equilibrate(ChemicalState& state, const Partition& partition, const EquilibriumOptions& options)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    # find a state based on the equilibrium setup
    state = equilibrate(problem)

    # change pressure and temperature
    state.setTemperature(problem.temperature() + 5)
    state.setPressure(problem.pressure() + 5)

    options = EquilibriumOptions()
    partition = Partition(system)

    # compute equilibrium for state with new temperature and pressure
    equilibriumResult = equilibrate(state, partition, options)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume"
            )
        ),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite at 60 °C and 300 bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved at 60 °C and 300 bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrate_overload_7(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of a problem using
    equilibrate(ChemicalState& state, const EquilibriumProblem& problem)
    and
    equilibrate(ChemicalState& state, const EquilibriumInverseProblem& problem)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    state = ChemicalState(system)

    equilibriumResult = equilibrate(state, problem)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts"
            )
        ),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume"
            )
        ),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite at 60 °C and 300 bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved at 60 °C and 300 bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrate_overload_8(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of a problem using
    equilibrate(ChemicalState& state, const EquilibriumProblem& problem, const EquilibriumOptions& options)
    and
    equilibrate(ChemicalState& state, const EquilibriumInverseProblem& problem, const EquilibriumOptions& options)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    state = ChemicalState(system)

    options = EquilibriumOptions()

    equilibriumResult = equilibrate(state, problem, options)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14))
