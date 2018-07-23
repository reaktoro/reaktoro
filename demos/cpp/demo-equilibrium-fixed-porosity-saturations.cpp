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

int main()
{
    ChemicalEditor editor;
    editor.addAqueousPhase("H O Na Cl C Ca Mg Si");
    editor.addGaseousPhase({"H2O(g)", "CO2(g)"});
    editor.addMineralPhase("Calcite");
    editor.addMineralPhase("Magnesite");
    editor.addMineralPhase("Dolomite");
    editor.addMineralPhase("Quartz");

    ChemicalSystem system(editor);

    EquilibriumOptions options;
    options.hessian = GibbsHessian::Exact;

    EquilibriumCompositionProblem composition(system);
    composition.setTemperature(100, "celsius");
    composition.setPressure(60, "bar");
    composition.setAqueousComposition("1 molal NaCl");
    composition.setGaseousComposition("0.90 CO2; 0.10 H2O");
    composition.setSolidComposition("0.10 Calcite; 0.05 Magnesite; 0.05 Dolomite; 0.80 Quartz");
    composition.setAqueousSaturation(0.80);
    composition.setGaseousSaturation(0.20);
    composition.setPorosity(0.3);

    ChemicalState state = equilibrate(composition, options);

    std::cout << state << std::endl;
}
