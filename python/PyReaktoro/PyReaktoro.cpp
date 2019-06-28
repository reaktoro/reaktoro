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

#include <PyReaktoro/PyReaktoro.hpp>
using namespace Reaktoro;

PYBIND11_MODULE(PyReaktoro, m)
{
    // Common module
    exportAutoDiff(m);
    exportIndex(m);
    exportMatrix(m);
    exportOutputter(m);
    exportReactionEquation(m);
    exportStringList(m);
    exportUnits(m);

    // Core module
    exportChemicalOutput(m);
    exportChemicalPlot(m);
    exportChemicalProperties(m);
    exportChemicalProperty(m);
    exportChemicalQuantity(m);
    exportChemicalState(m);
    exportChemicalSystem(m);
    exportConnectivity(m);
    exportElement(m);
    exportPartition(m);
    exportPhase(m);
    exportReaction(m);
    exportReactionSystem(m);
    exportSpecies(m);
    exportThermoProperties(m);
    exportUtils(m);

    // Equilibrium module
    exportEquilibriumCompositionProblem(m);
    exportEquilibriumInverseProblem(m);
    exportEquilibriumOptions(m);
    exportEquilibriumPath(m);
    exportEquilibriumProblem(m);
    exportEquilibriumResult(m);
    exportEquilibriumSensitivity(m);
    exportEquilibriumSolver(m);
    exportEquilibriumUtils(m);
    exportSmartEquilibriumOptions(m);
    exportSmartEquilibriumResult(m);
    exportSmartEquilibriumSolver(m);

    // Backends module
    exportInterface(m); // *** Warning *** exportInterface must be called before exportGems, exportPhreeqc, etc.
    exportGems(m);
    exportPhreeqc(m);
    exportPhreeqcEditor(m);

    // Interpreter module
    exportInterpreter(m);

    // Kinetics module
    exportKineticOptions(m);
    exportKineticPath(m);
    exportKineticSolver(m);

    // Math module
    exportODE(m);

    // Optimization module
    exportNonlinearOptions(m);
    exportOptimumMethod(m);
    exportOptimumOptions(m);
    exportOptimumResult(m);
    exportOptimumState(m);

    // Reactions module
    exportMineralCatalyst(m);
    exportMineralMechanism(m);
    exportMineralReaction(m);

    // Thermodynamics module
    exportStateOfMatter(m);
    exportDatabase(m);
    exportChemicalEditor(m);
    exportThermo(m);
    exportAqueousChemicalModelDebyeHuckel(m);
    exportAqueousPhase(m);
    exportGaseousPhase(m);
    exportMineralPhase(m);
    exportAqueousSpecies(m);
    exportGaseousSpecies(m);
    exportMineralSpecies(m);
    exportWater(m);

    // Transport module
    exportChemicalField(m);
    exportMesh(m);
    exportReactiveTransportOptions(m);
    exportReactiveTransportResult(m);
    exportReactiveTransportSolver(m);
    exportTransportSolver(m);
}
