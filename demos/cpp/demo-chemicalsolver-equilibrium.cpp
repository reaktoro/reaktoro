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
//
//#include <Reaktoro/Reaktoro.hpp>
//using namespace Reaktoro;
//
//int main()
//{
//    Index npoints = 10;
//
//    ChemicalEditor editor;
//    editor.addAqueousPhaseWithElements("H O Na Cl C Ca Mg");
//    editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
//    editor.addMineralPhase("Calcite");
//    editor.addMineralPhase("Dolomite");
//
//    ChemicalSystem system(editor);
//
//    EquilibriumCompositionProblem composition(system);
//    composition.setAqueousComposition("1 molal NaCl");
//    composition.setGaseousComposition("CO2");
//    composition.setSolidComposition("0.1 Calcite; 0.9 Dolomite");
//    composition.setAqueousSaturation(0.8);
//    composition.setGaseousSaturation(0.2);
//    composition.setPorosity(0.3);
//
//    ChemicalState state = equilibrate(composition);
//
//    ChemicalSolver solver(system, npoints);
//
//    Index Ee = system.numElements();
//
//    Vector T = constants(npoints, state.temperature());
//    Vector P = constants(npoints, state.pressure());
//    Matrix be(Ee, npoints);
//    be.colwise() = state.elementAmounts();
//
//    solver.equilibrate(T, P, be);
//
//    for(auto& state : solver.states())
//        std::cout << state << std::endl;
//
//    std::cout << "porosity = \n" << solver.porosity() << std::endl;
//    std::cout << "densities[0] = \n" << solver.fluidDensities()[0] << std::endl;
//    std::cout << "densities[1] = \n" << solver.fluidDensities()[1] << std::endl;
//    std::cout << "saturations[0] = \n" << solver.fluidSaturations()[0] << std::endl;
//    std::cout << "saturations[1] = \n" << solver.fluidSaturations()[1] << std::endl;
//}

int main(int argc, char **argv) {
    // TODO implement the above demo using the new design of Reaktoro
}
