// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

const auto NUMBER_OF_STEPS = 50;
const auto WARMSTART = true;

auto create_system() -> ChemicalSystem
{
    auto database = Database("supcrt07.xml");

    auto editor = ChemicalEditor(database);
    editor.addAqueousPhaseWithElementsOf({"C", "Ca", "Cl", "Fe", "H", "Na", "O", "S", "Ba", "Sr"});
    editor.addMineralPhase("Anhydrite");
    editor.addMineralPhase("Barite");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Celestite");
    editor.addMineralPhase("Siderite");
    editor.addMineralPhase("Pyrrhotite");

    return ChemicalSystem(editor);
}

auto create_default_problem(const ChemicalSystem& system) -> EquilibriumProblem
{
    const auto pressure = 13660337.317389656; // Pa
    const auto temperature = 64;              // C

    auto problem = EquilibriumProblem(system);
    problem.setTemperature(temperature, "celsius");
    problem.setPressure(pressure, "Pa");
    problem.add("H2O(l)", 394.76384905586252, "kg");
    problem.add("CO2(aq)", 25.350676739651345, "mol");
    problem.add("HCO3-", 40.7007440302866, "mol");
    problem.add("CO3--", 0.026502456032096638, "mol");
    problem.add("H2S(aq)", 5.7146433306057664e-24, "mol");
    problem.add("HS-", 6.0271635093622933e-24, "mol");
    problem.add("H+", 0.00026713431941822389, "mol");
    problem.add("OH-", 0.00070239234495091823, "mol");
    problem.add("SO4--", 0.006002566648441959, "kg");
    problem.add("Ca++", 0.011204791077091658, "kg");
    problem.add("Na+", 3.1086937041933971, "kg");
    problem.add("Cl-", 3.3645158450242447, "kg");
    problem.add("Fe++", 0.00023317970573647533, "kg");
    problem.add("Ba++", 0.0080034221979226137, "kg");
    problem.add("Sr++", 0.017607528835429747, "kg");
    problem.add("Calcite", 100000000, "kg");
    problem.add("Siderite", 0.01795172676219817, "mol");

    return problem;
}

// Using only approximate, the problem doens"t happen.
auto test_with_same_solver_instance() -> void
{
    auto system = create_system();
    auto default_problem = create_default_problem(system);
    auto state = ChemicalState(system);
    auto solver = EquilibriumSolver(system);
    auto options = EquilibriumOptions();
    options.warmstart = WARMSTART;
    options.optimum.tolerance = 1e-10;
    solver.setOptions(options);

    for(auto i = 0; i < NUMBER_OF_STEPS; ++i) {
        std::cout << "test_with_same_solver_instance:" << i << "/" << NUMBER_OF_STEPS << std::endl;
        solver.solve(state, default_problem);
        default_problem.addSpecies("H2S(aq)", 0.1, "mg");
    }

    state.output("UsingSameSolverInstance.txt");
}

auto test_with_different_solver_instance() -> void
{
    auto system = create_system();
    auto default_problem = create_default_problem(system);
    auto state = ChemicalState(system);
    auto options = EquilibriumOptions();
    options.warmstart = WARMSTART;
    options.optimum.tolerance = 1e-10;

    for(auto i = 0; i < NUMBER_OF_STEPS; ++i) {
        std::cout << "test_with_different_solver_instance:" << i << "/" << NUMBER_OF_STEPS << std::endl;
        auto solver = EquilibriumSolver(system);
        solver.setOptions(options);

        solver.solve(state, default_problem);
        default_problem.addSpecies("H2S(aq)", 0.1, "mg");
    }

    state.output("UsingDifferentSolverInstance.txt");
}

int main()
{
    test_with_same_solver_instance();
    test_with_different_solver_instance();
}
