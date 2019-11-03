# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2018 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.

from reaktoro import EquilibriumSolver, ChemicalState, EquilibriumProblem, equilibrate


def _create_equilibrium_problem(partition_with_inert_gaseous_phase):
    # Equilibrate 1 kg of H2O and 1 mol of CO2, but CO2(g) and H2O(g) are inert
    problem = EquilibriumProblem(partition_with_inert_gaseous_phase)
    problem.add("H2O", 1, "kg")
    problem.add("CO2", 1, "mol")

    return problem


def _create_chemical_state(chemical_system):
    # Set the initial amounts of CO2(g) and H2O(g)
    state = ChemicalState(chemical_system)

    # Set the initial amounts of CO2(g) and H2O(g)
    state.setSpeciesAmount("CO2(g)", 1, "mol")
    state.setSpeciesAmount("H2O(g)", 1, "mmol")

    return state


def test_equilibrium_solver_with_equilibrate(partition_with_inert_gaseous_phase, chemical_system):
    problem = _create_equilibrium_problem(partition_with_inert_gaseous_phase)
    state = _create_chemical_state(chemical_system)

    # Compute equilibrium state; the amounts of CO2(g) and H2O(g) should remain the same
    equilibrate(state, problem)

    # Assert the amounts of CO2(g) and H2O(g) are the same as initially set
    assert state.speciesAmount("CO2(g)") == 1.0
    assert state.speciesAmount("H2O(g)") == 0.001


def test_equilibrium_solver_with_equilibrium_problem(
    partition_with_inert_gaseous_phase, chemical_system
):
    problem = _create_equilibrium_problem(partition_with_inert_gaseous_phase)
    state = _create_chemical_state(chemical_system)

    # Compute equilibrium state; the amounts of CO2(g) and H2O(g) should remain the same
    solver = EquilibriumSolver(chemical_system)
    solver.setPartition(problem.partition())
    solver.solve(state, problem)

    # Assert the amounts of CO2(g) and H2O(g) are the same as initially set
    assert state.speciesAmount("CO2(g)") == 1.0
    assert state.speciesAmount("H2O(g)") == 0.001


def test_equilibrium_solver_with_equilibrium_problem_adding_argon(
    partition_with_inert_gaseous_phase_adding_argon, chemical_system_adding_argon
):
    problem = _create_equilibrium_problem(partition_with_inert_gaseous_phase_adding_argon)
    state = _create_chemical_state(chemical_system_adding_argon)

    # Compute equilibrium state; the amounts of CO2(g) and H2O(g) should remain the same
    solver = EquilibriumSolver(chemical_system_adding_argon)
    solver.setPartition(problem.partition())
    solver.solve(state, problem)

    # Assert the amounts of CO2(g) and H2O(g) are the same as initially set
    assert state.speciesAmount("CO2(g)") == 1.0
    assert state.speciesAmount("H2O(g)") == 0.001


def test_equilibrium_solver_with_chemical_state(
    partition_with_inert_gaseous_phase, chemical_system
):
    state = _create_chemical_state(chemical_system)

    # Compute equilibrium state; the amounts of CO2(g) and H2O(g) should remain the same
    solver = EquilibriumSolver(chemical_system)
    solver.setPartition(partition_with_inert_gaseous_phase)

    # Set the amounts of H2O(l) and CO2(aq)
    state.setSpeciesMass("H2O(l)", 55, "kg")
    state.setSpeciesAmount("CO2(aq)", 1, "mol")

    # Check inert elements are taken care of when solving with state input
    solver.solve(state)

    # Assert the amounts of CO2(g) and H2O(g) are the same as initially set
    assert state.speciesAmount("CO2(g)") == 1.0
    assert state.speciesAmount("H2O(g)") == 0.001
