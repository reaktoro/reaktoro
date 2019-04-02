import numpy as np
import pytest

from reaktoro import ChemicalState, equilibrate, EquilibriumSolver, EquilibriumOptions


@pytest.mark.parametrize(
    "setup",
    [
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar"
            )
        ),
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
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass"
            )
        ),
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
        "Eq Prob-H2O CO2 NaCl and Halite 60C 300bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrium_solver_solve_overload_1(setup, state_regression):
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

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize(
    "setup",
    [
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar"
            )
        ),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite 60C 300bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
    ],
)
def test_equilibrium_solver_solve_overload_2(setup, state_regression):
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

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize(
    "setup",
    [
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar"
            )
        ),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite 60C 300bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
    ],
)
def test_equilibrium_solver_solve_overload_3(setup, state_regression):
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

    solver.solve(
        state, problem.temperature(), problem.pressure(), problem.elementAmounts()
    )

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize(
    "setup",
    [
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar"
            )
        ),
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
        (
            pytest.lazy_fixture(
                "equilibrium_inverse_with_h2o_nacl_caco3_calcilte_and_fixed_mass"
            )
        ),
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
        "Eq Prob-H2O CO2 NaCl and Halite 60C 300bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
        "InvEq Prob-H O Na Cl Ca Mg C with fixed amount and activity",
        "InvEq Prob-H O Na Cl Ca Mg C with defined pH",
        "InvEq Prob-H O Na Cl Ca C Calcite with defined pH and fixed amount",
        "InvEq Prob-H2O CO2 NaCl CaCO3 Calcite with fixed species mass",
        "InvEq Prob-H2O NaCl CaCO3 CO2 with fixed mass amount and alkalinity",
        "InvEq Prob-H2H NaCl CaCO3 CO2 Calcite with fixed phase volume",
    ],
)
def test_equilibrium_solver_approx_overload_1(setup, state_regression):
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
    state.setTemperature(problem.temperature() + 15)
    state.setPressure(problem.pressure() + 20)

    solver = EquilibriumSolver(system)

    solver.approximate(state)

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-16))

@pytest.mark.parametrize(
    "setup",
    [
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar"
            )
        ),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite 60C 300bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
    ],
)
def test_equilibrium_solver_approx_overload_2(setup, state_regression):
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

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-16))


@pytest.mark.parametrize(
    "setup",
    [
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_60C_300bar"
            )
        ),
        (
            pytest.lazy_fixture(
                "equilibrium_problem_with_h2o_co2_nacl_halite_dissolved_60C_300bar"
            )
        ),
        (pytest.lazy_fixture("equilibrium_problem_with_h2o_feoh2_feoh3_nh3_magnetite")),
    ],
    ids=[
        "Eq Prob-H2O CO2 NaCl and Halite 60C 300bar",
        "Eq Prob-H2O CO2 NaCl and Halite already dissolved 60C 300bar",
        "Eq Prob-H2O Fe(OH)2 Fe(OH)3 NH3 and Magnetite",
    ],
)
def test_equilibrium_solver_approx_overload_3(setup, state_regression):
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

    solver.approximate(
        state, problem.temperature(), problem.pressure(), problem.elementAmounts()
    )

    state_regression.check(state, default_tol=dict(atol=1e-5, rtol=1e-16))
