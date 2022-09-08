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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

void exportActivityModel(py::module& m);
void exportActivityProps(py::module& m);
void exportAggregateState(py::module& m);
void exportChemicalFormula(py::module& m);
void exportChemicalProps(py::module& m);
void exportChemicalPropsPhase(py::module& m);
void exportChemicalState(py::module& m);
void exportChemicalSystem(py::module& m);
void exportData(py::module& m);
void exportDatabase(py::module& m);
void exportElement(py::module& m);
void exportElementalComposition(py::module& m);
void exportElementList(py::module& m);
void exportEmbedded(py::module& m);
void exportFormationReaction(py::module& m);
void exportPhase(py::module& m);
void exportPhaseList(py::module& m);
void exportPhases(py::module& m);
void exportParam(py::module& m);
void exportParams(py::module& m);
void exportReaction(py::module& m);
void exportReactions(py::module& m);
void exportReactionEquation(py::module& m);
void exportReactionList(py::module& m);
void exportReactionThermoProps(py::module& m);
void exportReactionStandardThermoModel(py::module& m);
void exportReactionStandardThermoProps(py::module& m);
void exportSpecies(py::module& m);
void exportSpeciesList(py::module& m);
void exportSpeciesThermoProps(py::module& m);
void exportStandardThermoModel(py::module& m);
void exportStandardThermoProps(py::module& m);
void exportStateOfMatter(py::module& m);
void exportThermoProps(py::module& m);
void exportThermoPropsPhase(py::module& m);
void exportCoreUtils(py::module& m);

void exportCore(py::module& m)
{
    exportActivityModel(m);
    exportActivityProps(m);
    exportAggregateState(m);
    exportChemicalFormula(m);
    exportChemicalProps(m);
    exportChemicalPropsPhase(m);
    exportChemicalState(m);
    exportChemicalSystem(m);
    exportData(m);
    exportDatabase(m);
    exportElement(m);
    exportElementalComposition(m);
    exportElementList(m);
    exportEmbedded(m);
    exportFormationReaction(m);
    exportPhase(m);
    exportPhaseList(m);
    exportPhases(m);
    exportParam(m);
    exportParams(m);
    exportReaction(m);
    exportReactions(m);
    exportReactionEquation(m);
    exportReactionList(m);
    exportReactionThermoProps(m);
    exportReactionStandardThermoProps(m);
    exportReactionStandardThermoModel(m);
    exportSpecies(m);
    exportSpeciesList(m);
    exportSpeciesThermoProps(m);
    exportStandardThermoProps(m);
    exportStandardThermoModel(m);
    exportStateOfMatter(m);
    exportThermoProps(m);
    exportThermoPropsPhase(m);
    exportCoreUtils(m);

}
