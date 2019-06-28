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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
namespace py = pybind11;

namespace Reaktoro {

// Common module
void exportAutoDiff(py::module& m);
void exportEigen(py::module& m);
void exportIndex(py::module& m);
void exportMatrix(py::module& m);
void exportOutputter(py::module& m);
void exportReactionEquation(py::module& m);
void exportStandardTypes(py::module& m);
void exportStringList(py::module& m);
void exportUnits(py::module& m);

// Core module
void exportChemicalOutput(py::module& m);
void exportChemicalPlot(py::module& m);
void exportChemicalProperties(py::module& m);
void exportChemicalProperty(py::module& m);
void exportChemicalQuantity(py::module& m);
void exportChemicalState(py::module& m);
void exportChemicalSystem(py::module& m);
void exportConnectivity(py::module& m);
void exportElement(py::module& m);
void exportPartition(py::module& m);
void exportPhase(py::module& m);
void exportReaction(py::module& m);
void exportReactionSystem(py::module& m);
void exportSpecies(py::module& m);
void exportThermoProperties(py::module& m);
void exportUtils(py::module& m);

// Equilibrium module
void exportEquilibriumCompositionProblem(py::module& m);
void exportEquilibriumInverseProblem(py::module& m);
void exportEquilibriumOptions(py::module& m);
void exportEquilibriumPath(py::module& m);
void exportEquilibriumProblem(py::module& m);
void exportEquilibriumResult(py::module& m);
void exportEquilibriumSensitivity(py::module& m);
void exportEquilibriumSolver(py::module& m);
void exportEquilibriumUtils(py::module& m);
void exportSmartEquilibriumSolver(py::module& m);

// Backends module
void exportGems(py::module& m);
void exportInterface(py::module& m);
void exportPhreeqc(py::module& m);
void exportPhreeqcEditor(py::module& m);

// Interpreter module
void exportInterpreter(py::module& m);

// Kinetics module
void exportKineticOptions(py::module& m);
void exportKineticPath(py::module& m);
void exportKineticSolver(py::module& m);

// Math module
void exportODE(py::module& m);

// Optimization module
void exportNonlinearOptions(py::module& m);
void exportOptimumMethod(py::module& m);
void exportOptimumOptions(py::module& m);
void exportOptimumResult(py::module& m);
void exportOptimumState(py::module& m);

// Reactions module
void exportMineralCatalyst(py::module& m);
void exportMineralMechanism(py::module& m);
void exportMineralReaction(py::module& m);

// Thermodynamics module
void exportStateOfMatter(py::module& m);
void exportChemicalEditor(py::module& m);
void exportDatabase(py::module& m);
void exportThermo(py::module& m);
void exportAqueousChemicalModelDebyeHuckel(py::module& m);
void exportAqueousPhase(py::module& m);
void exportGaseousPhase(py::module& m);
void exportLiquidPhase(py::module& m);
void exportMineralPhase(py::module& m);
void exportAqueousSpecies(py::module& m);
void exportGaseousSpecies(py::module& m);
void exportLiquidSpecies(py::module& m);
void exportMineralSpecies(py::module& m);
void exportWater(py::module& m);
void exportOil(py::module& m);

// Transport module
void exportChemicalField(py::module& m);
void exportMesh(py::module& m);
void exportTransportSolver(py::module& m);
void exportReactiveTransportSolver(py::module& m);

} // namespace Reaktoro
