// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

int main()
{
    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution("H2O(aq) H+ OH- Ca+2 HCO3- CO3-2 CO2(aq)");
    solution.set(AcMoDav);

    MineralPhase calcite("Calcite");

    Phases phases(db);
    phases.add(solution);
    phases.add(calcite);

    Params params = Params::embedded("PalandriKharaka.yaml");

    MineralReactions reactions("Calcite");
    reactions.setRateModel(ReactionRateModelPalandriKharaka(params));

    ChemicalSystem system(phases, reactions);

    ChemicalState state(system);
    state.set("H2O(aq)", 1.0, "kg");
    state.set("Calcite", 1.0, "mol");
    state.surfaceArea("Calcite", 1.0, "m2");

    KineticsSolver solver(system);

    real dt = 0.1;

    EquilibriumOptions options;
    options.optima.output.active = true;

    solver.setOptions(options);

    EquilibriumResult result;

    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);
    result = solver.solve(state, dt);

    std::cout << state << std::endl;

    ChemicalState state_eq(system);
    state_eq.set("H2O(aq)", 1.0, "kg");
    state_eq.set("Calcite", 1.0, "mol");
    equilibrate(state_eq);

    std::cout << state_eq << std::endl;
}
