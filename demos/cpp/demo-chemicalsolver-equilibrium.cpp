// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    Index npoints = 10;

    ChemicalEditor editor;
    editor.addAqueousPhase("H O Na Cl C Ca Mg");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Dolomite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.add("H2O", 1.0, "kg");
    problem.add("NaCl", 1.0, "mol");
    problem.add("CaCO3", 2.0, "mol");
    problem.add("MgCO3", 1.0, "mol");

    ChemicalState state = equilibrate(problem);

    std::cout << "V(fluid) = " << state.fluidVolume().val << std::endl;
    std::cout << "V(solid) = " << state.solidVolume().val << std::endl;

    return 0;

    ChemicalSolver solver(system, npoints);

    Partition partition(system);
    partition.setFluidPhases({"Aqueous"});
    solver.setPartition(partition);

    Index Ee = partition.numEquilibriumElements();

    Vector T = constants(npoints, 300.0);  // 300 K
    Vector P = constants(npoints, 60.0e5); //  60 bar
    Matrix be(Ee, npoints);
    be.colwise() = problem.elementAmounts();

    solver.equilibrate(T.data(), P.data(), be.data());

    for(auto& state : solver.states())
        std::cout << state << std::endl;

    ChemicalField porosity;

    solver.porosity(porosity);

    for(Index i = 0; i < npoints; ++i)
        std::cout << porosity.val[i] << std::endl;
}
