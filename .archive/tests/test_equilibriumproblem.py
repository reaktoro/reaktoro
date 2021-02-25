# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright (C) 2014-2021 Allan Leal
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

from reaktoro import ChemicalState, EquilibriumProblem


def test_equilibrium_problem_add_by_chemical_state(partition_with_inert_gaseous_phase, chemical_system):
    state = ChemicalState(chemical_system)
    state.setSpeciesAmount("CO2(g)", 10.0)

    problem = EquilibriumProblem(partition_with_inert_gaseous_phase)
    problem.add(state)

    for element in problem.elementAmounts():
        assert element == 0.0
    assert problem.partition().numInertSpecies() == 2


def test_equilibrium_problem_add(partition_with_inert_gaseous_phase):
    problem = EquilibriumProblem(partition_with_inert_gaseous_phase)
    problem.add("CO2", 10.0, 'mol')

    assert sum(problem.elementAmounts()) == 30.0
    assert problem.partition().numInertSpecies() == 2
