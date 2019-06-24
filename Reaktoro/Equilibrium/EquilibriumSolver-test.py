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

import pytest

from reaktoro import EquilibriumSolver, ChemicalState, EquilibriumProblem


@pytest.mark.xfail(reason='fix equilibrium solver state with only inert species')
def test_equilibrium_solver_with_inert_gaseous_phase_using_chemical_state(partition_with_inert_gaseous_phase, chemical_system):
    solver = EquilibriumSolver(partition_with_inert_gaseous_phase)

    state = ChemicalState(chemical_system)
    state.setSpeciesAmount("CO2(g)", 10.0)
    element_amounts_before_solver = state.elementAmounts()
    solver.approximate(state)
    for (element_before, element_after) in zip(element_amounts_before_solver, state.elementAmounts()):
        assert element_before == pytest.approx(element_after)


@pytest.mark.xfail(reason='fix equilibrium solver with inert species')
def test_equilibrium_solver_with_inert_gaseous_phase_using_equilibrium_problem(partition_with_inert_gaseous_phase, chemical_system):
    problem = EquilibriumProblem(partition_with_inert_gaseous_phase)
    solver = EquilibriumSolver(partition_with_inert_gaseous_phase)

    state = ChemicalState(chemical_system)

    state.setSpeciesAmount("CO2(g)", 10.0)
    state.setSpeciesAmount("CO2(aq)", 10.0)
    problem.add(state)

    element_amounts_before_solver = problem.elementAmounts()
    solver.approximate(state, problem)

    for (element_before, element_after) in zip(element_amounts_before_solver, state.elementAmounts()):
        assert element_before == pytest.approx(element_after)
