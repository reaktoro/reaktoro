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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Allan Leal (16 July 2021)
//
// and since revised by:
//   • Allan Leal (22 July 2021)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution(speciate("Na Cl C Ca Mg Si"), exclude("organic"));
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    GaseousPhase gases("CO2(g) H2O(g)");
    gases.setActivityModel(ActivityModelPengRobinson());

    MineralPhases minerals("Halite Calcite Magnesite Dolomite Quartz");

    ChemicalSystem system(db, solution, gases, minerals);

    ChemicalState statex(system);
    statex.temperature(60.0, "celsius");
    statex.pressure(100.0, "bar");
    statex.set("H2O(aq)"  , 1.00, "kg");
    statex.set("Halite"   , 1.00, "mol");
    statex.set("Calcite"  , 1.00, "mol");
    statex.set("Magnesite", 1.00, "mol");
    statex.set("Quartz"   , 1.00, "mol");

    equilibrate(statex);

    ChemicalProps propsx(statex);
    const auto Vx = propsx.volume();
    const auto Ux = propsx.internalEnergy();

    statex.output("state-expected.txt");
    propsx.output("props-expected.txt");

    EquilibriumSpecs specs(system);

    auto idxV = specs.addInput("V");
    auto idxU = specs.addInput("U");

    ConstraintEquation volumeConstraint;
    volumeConstraint.id = "VolumeConstraint";
    volumeConstraint.fn = [=](const ChemicalState& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().volume() - w[idxV];
    };

    ConstraintEquation internalEnergyConstraint;
    internalEnergyConstraint.id = "InternalEnergyConstraint";
    internalEnergyConstraint.fn = [=](const ChemicalState& state, VectorXrConstRef p, VectorXrConstRef w)
    {
        return state.props().internalEnergy() - w[idxU];
    };

    specs.addConstraint(volumeConstraint);
    specs.addConstraint(internalEnergyConstraint);

    EquilibriumConditions conditions(specs);
    conditions.set("V", Vx);
    conditions.set("U", Ux);
    conditions.setLowerBoundPressure(1.0, "bar");

    ChemicalState state(system);
    state.temperature(25.0, "celsius");
    state.pressure(1.0, "bar");
    state.set("H2O(aq)"  , 1.00, "kg");
    state.set("Halite"   , 1.00, "mol");
    state.set("Calcite"  , 1.00, "mol");
    state.set("Magnesite", 1.00, "mol");
    state.set("Quartz"   , 1.00, "mol");

    EquilibriumSolver solver(specs);

    solver.solve(state, conditions);

    ChemicalProps props(state);

    state.output("state.txt");
    props.output("props.txt");

    // Check if props.txt and props-expected.txt are numerically equivalent.
    std::cout << "Success! Check outputted files `state.txt`, `props.txt`, `state-expected.txt`, `props-expected.txt`." << std::endl;
    std::cout << "Verify if `props.txt` and `props-expected.txt` are numerically equivalent." << std::endl;

    return 0;
}
