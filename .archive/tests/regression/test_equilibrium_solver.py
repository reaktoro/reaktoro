import numpy as np
import pytest

from reaktoro import ChemicalState, equilibrate, EquilibriumSolver, EquilibriumOptions


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (pytest.lazy_fixture("equilibrium_problem_using_thermofun_aq17_database")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume")),
    ],
    ids=[
        "H2O CO2 NaCl and Halite 60C 300bar",
        "H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "ThermoFun database aq17",
        "InvEq-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq-H O Na Cl Ca Mg C with defined pH",
        "InvEq-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_solve_overload_1(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of an equilibrium of a state using
    EquilibriumSolver::solve(ChemicalState& state)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    state = equilibrate(problem)
    state.setTemperature(problem.temperature() + 30)
    state.setPressure(problem.pressure() + 40)

    solver = EquilibriumSolver(system)

    solver.solve(state)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (pytest.lazy_fixture("equilibrium_problem_using_thermofun_aq17_database")),
    ],
    ids=[
        "H2O CO2 NaCl and Halite 60C 300bar",
        "H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "ThermoFun database aq17",
    ],
)
def test_solve_overload_2(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of an equilibrium of a state using
    EquilibriumSolver::solve(ChemicalState& state, const EquilibriumProblem& problem)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.solve(state, problem)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (pytest.lazy_fixture("equilibrium_problem_using_thermofun_aq17_database")),
    ],
    ids=[
        "H2O CO2 NaCl and Halite 60C 300bar",
        "H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "ThermoFun database aq17",
    ],
)
def test_solve_overload_3(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of an equilibrium of a state using
    EquilibriumSolver::solve(ChemicalState& state, double T, double P, VectorConstRef be)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    T = problem.temperature()
    P = problem.pressure()
    b = problem.elementAmounts()

    solver.solve(state, T, P, b)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14))


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (pytest.lazy_fixture("equilibrium_problem_using_thermofun_aq17_database")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_fixed_amount_and_activity")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_mg_c_defined_ph")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h_o_na_cl_ca_c_calcite_ph_and_fixed_amounts")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_co2_fixed_mass_amount_and_alkalinity")),
        (pytest.lazy_fixture("equilibrium_inverse_with_h2o_nacl_caco3_co2_calcite_fixed_phase_volume")),
    ],
    ids=[
        "H2O CO2 NaCl and Halite 60C 300bar",
        "H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "ThermoFun database aq17",
        "InvEq-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq-H O Na Cl Ca Mg C with defined pH",
        "InvEq-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_approx_overload_1(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of an equilibrium of a state using
    EquilibriumSolver::approximate(ChemicalState& state)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    state = equilibrate(problem)
    state.setTemperature(problem.temperature() + 10)
    state.setPressure(problem.pressure() + 10)

    solver = EquilibriumSolver(system)

    solver.approximate(state)

    exclude = ['pH [-]']

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14), exclude=exclude)

@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (pytest.lazy_fixture("equilibrium_problem_using_thermofun_aq17_database")),
    ],
    ids=[
        "H2O CO2 NaCl and Halite 60C 300bar",
        "H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "ThermoFun database aq17",
    ],
)
def test_approx_overload_2(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of an equilibrium of a state using
    EquilibriumSolver::approximate(ChemicalState& state, const EquilibriumProblem& problem)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    solver.approximate(state, problem)

    exclude = ['pH [-]']

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14), exclude=exclude)


@pytest.mark.parametrize(
    "setup",
    [
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar")),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
        (pytest.lazy_fixture("equilibrium_problem_using_thermofun_aq17_database")),
    ],
    ids=[
        "H2O CO2 NaCl and Halite 60C 300bar",
        "H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "ThermoFun database aq17",
    ],
)
def test_approx_overload_3(setup, state_regression):
    """
    An integration test that checks result's reproducibility of
    the calculation of an equilibrium of a state using
    EquilibriumSolver::approximate(ChemicalState& state, double T, double P, VectorConstRef be)
    @param setup
        a tuple that has some objects from problem setup
        (system, problem)
    """
    (system, problem) = setup

    state = ChemicalState(system)

    solver = EquilibriumSolver(system)

    T = problem.temperature()
    P = problem.pressure()
    b = problem.elementAmounts()

    solver.approximate(state, T, P, b)

    exclude = ['pH [-]']

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-14), exclude=exclude)
