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
    editor.addGaseousPhase("H2O(g) CO2(g)");
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Dolomite");

    ChemicalSystem system(editor);

    EquilibriumProblem problem(system);
    problem.add("H2O", 1.0, "kg");
    problem.add("NaCl", 1.0, "mol");
    problem.add("CaCO3", 2.0, "mol");
    problem.add("MgCO3", 1.0, "mol");
    problem.add("CO2", 3.0, "mol");

    ChemicalState state = equilibrate(problem);

    state.scalePhaseVolume("Aqueous", 0.4);
    state.scalePhaseVolume("Gaseous", 0.1);
    state.scaleSolidVolume(0.5);

    std::cout << state << std::endl;

    ChemicalProperties properties = state.properties();
    std::cout << "V(fluid) = " << properties.fluidVolume().val << std::endl;
    std::cout << "V(solid) = " << properties.solidVolume().val << std::endl;
//    std::cout << "phi = " << state.porosity().val << std::endl;
//    std::cout << "sat = \n" << state.saturations().val << std::endl;
    std::cout << "rho = \n" << properties.phaseDensities().val << std::endl;

    ChemicalSolver solver(system, npoints);

    Index Ee = system.numElements();

    Vector T = constants(npoints, state.temperature());
    Vector P = constants(npoints, state.pressure());
    Matrix be(Ee, npoints);
    be.colwise() = state.elementAmounts();

    solver.equilibrate(T.data(), P.data(), be.data());

    for(auto& state : solver.states())
        std::cout << state << std::endl;


    auto porosity = solver.porosity();

    std::cout << "porosity = \n" << porosity.val() << std::endl;
}
