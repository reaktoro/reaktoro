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

#include <pybind11/stl_bind.h>

namespace Reaktoro {

// Common module
extern void exportAutoDiff(py::module& m);
extern void exportEigen(py::module& m);
extern void exportIndex(py::module& m);
extern void exportOpenlibm(py::module& m);
extern void exportMatrix(py::module& m);
extern void exportOutputter(py::module& m);
extern void exportReactionEquation(py::module& m);
extern void exportStandardTypes(py::module& m);
extern void exportStringList(py::module& m);
extern void exportUnits(py::module& m);

// Core module
extern void exportChemicalOutput(py::module& m);
extern void exportChemicalPlot(py::module& m);
extern void exportChemicalProperties(py::module& m);
extern void exportChemicalProperty(py::module& m);
extern void exportChemicalQuantity(py::module& m);
extern void exportChemicalState(py::module& m);
extern void exportChemicalSystem(py::module& m);
extern void exportConnectivity(py::module& m);
extern void exportElement(py::module& m);
extern void exportPartition(py::module& m);
extern void exportPhase(py::module& m);
extern void exportReaction(py::module& m);
extern void exportReactionSystem(py::module& m);
extern void exportSpecies(py::module& m);
extern void exportThermoProperties(py::module& m);
extern void exportUtils(py::module& m);

// Equilibrium module
extern void exportEquilibriumCompositionProblem(py::module& m);
extern void exportEquilibriumInverseProblem(py::module& m);
extern void exportEquilibriumOptions(py::module& m);
extern void exportEquilibriumPath(py::module& m);
extern void exportEquilibriumProblem(py::module& m);
extern void exportEquilibriumResult(py::module& m);
extern void exportEquilibriumSensitivity(py::module& m);
extern void exportEquilibriumSolver(py::module& m);
extern void exportEquilibriumUtils(py::module& m);
extern void exportSmartEquilibriumSolver(py::module& m);

// Backends module
extern void exportGems(py::module& m);
extern void exportInterface(py::module& m);
extern void exportPhreeqc(py::module& m);
extern void exportPhreeqcEditor(py::module& m);

// Interpreter module
extern void exportInterpreter(py::module& m);

// Kinetics module
extern void exportKineticOptions(py::module& m);
extern void exportKineticPath(py::module& m);
extern void exportKineticSolver(py::module& m);

// Math module
extern void exportODE(py::module& m);
extern void exportBilinearInterpolator(py::module& m);

// Optimization module
extern void exportNonlinearOptions(py::module& m);
extern void exportOptimumMethod(py::module& m);
extern void exportOptimumOptions(py::module& m);
extern void exportOptimumResult(py::module& m);
extern void exportOptimumState(py::module& m);

// Reactions module
extern void exportMineralCatalyst(py::module& m);
extern void exportMineralMechanism(py::module& m);
extern void exportMineralReaction(py::module& m);

// Thermodynamics module
extern void exportStateOfMatter(py::module& m);
extern void exportChemicalEditor(py::module& m);
extern void exportDatabase(py::module& m);
extern void exportThermo(py::module& m);
extern void exportAqueousChemicalModelDebyeHuckel(py::module& m);
extern void exportAqueousPhase(py::module& m);
extern void exportFluidPhase(py::module& m);
extern void exportMineralPhase(py::module& m);
extern void exportAqueousSpecies(py::module& m);
extern void exportFluidSpecies(py::module& m);
extern void exportMineralSpecies(py::module& m);
extern void exportWater(py::module& m);
extern void exportThermoData(py::module& m);
extern void exportThermoDataProperties(py::module& m);

// Transport module
extern void exportChemicalField(py::module& m);
extern void exportMesh(py::module& m);
extern void exportTransportSolver(py::module& m);
extern void exportReactiveTransportSolver(py::module& m);

} // namespace Reaktoro

using namespace Reaktoro;

PYBIND11_MODULE(PyReaktoro, m)
{
    py::bind_vector<std::vector<double>>(m, "VectorDouble", "VectorDouble Descriptor");

    // Common module
    exportAutoDiff(m);
    exportIndex(m);
    exportOpenlibm(m);
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
    exportBilinearInterpolator(m);

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
    exportFluidPhase(m);
    exportMineralPhase(m);
    exportAqueousSpecies(m);
    exportFluidSpecies(m);
    exportMineralSpecies(m);
    exportWater(m);
    exportThermoData(m);
    exportThermoDataProperties(m);

    // Transport module
    exportChemicalField(m);
    exportMesh(m);
    exportTransportSolver(m);
    exportReactiveTransportSolver(m);
}
