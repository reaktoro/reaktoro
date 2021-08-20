// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Allan Leal (4 August 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

const auto T = 60.0 + 273.15; // temperature in K
const auto P = 10.0 * 1e5;    // pressure in Pa

auto computeSolubilityCaCO3(const ChemicalSystem& system) -> real
{
    ChemicalState state(system);
    state.setTemperature(T);
    state.setPressure(P);
    state.set("H2O(aq)", 1.0, "kg");
    state.set("Calcite", 1.0, "mol");

    equilibrate(state);

    AqueousProps aprops(state);

    return aprops.elementMolality("Ca");
}

int main()
{
    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution(speciate("C Ca"), exclude("organic"));
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    MineralPhase mineral("Calcite");

    ChemicalSystem system(db, solution, mineral);

    const auto solubilityCaCO3 = computeSolubilityCaCO3(system);

    ConstraintEquation constraint;
    constraint.id = "solubility[CaCO3]";
    constraint.fn = [=](const ChemicalProps& props, VectorXrConstRef w)
    {
        static AqueousProps aprops(system);
        aprops.update(props);
        return aprops.elementMolality("Ca") - solubilityCaCO3;
    };

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.addUnknownStandardChemicalPotential("Calcite");
    specs.addConstraint(constraint);

    ChemicalState state(system);
    state.setTemperature(T);
    state.setPressure(P);
    state.set("H2O(aq)", 1.0, "kg");
    state.set("Calcite", 10.0, "mol");

    EquilibriumSolver solver(specs);

    solver.solve(state);

    const auto G0_calcite_expected = system.species().get("Calcite").props(T, P).G0;
    const auto G0_calcite_computed = state.equilibrium().p()[0];

    std::cout << "=================================" << std::endl;
    std::cout << "G0(calcite) at 60 °C and 10 bar  " << std::endl;
    std::cout << "=================================" << std::endl;
    std::cout << "expected: " << G0_calcite_expected/1000.0 << " kJ/mol" << std::endl;
    std::cout << "computed: " << G0_calcite_computed/1000.0 << " kJ/mol" << std::endl;
    std::cout << "   error: " << abs((G0_calcite_computed - G0_calcite_expected)/G0_calcite_expected) * 100.0 << " %" << std::endl;
    std::cout << "=================================" << std::endl;

    return 0;
}
