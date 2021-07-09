// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

void exportReactionThermoModelConstLgK(py::module& m);
void exportReactionThermoModelGemsLgK(py::module& m);
void exportReactionThermoModelPhreeqcLgK(py::module& m);
void exportReactionThermoModelPressureCorrection(py::module& m);
void exportReactionThermoModelVantHoff(py::module& m);
void exportReactionThermoModelYAML(py::module& m);

void exportStandardThermoModelConstant(py::module& m);
void exportStandardThermoModelHKF(py::module& m);
void exportStandardThermoModelHollandPowell(py::module& m);
void exportStandardThermoModelMaierKelley(py::module& m);
void exportStandardThermoModelMineralHKF(py::module& m);
void exportStandardThermoModelWaterHKF(py::module& m);
void exportStandardThermoModelYAML(py::module& m);

void exportModels(py::module& m)
{
    exportReactionThermoModelConstLgK(m);
    exportReactionThermoModelGemsLgK(m);
    exportReactionThermoModelPhreeqcLgK(m);
    exportReactionThermoModelPressureCorrection(m);
    exportReactionThermoModelVantHoff(m);
    exportReactionThermoModelYAML(m);

    exportStandardThermoModelConstant(m);
    exportStandardThermoModelHKF(m);
    exportStandardThermoModelHollandPowell(m);
    exportStandardThermoModelMaierKelley(m);
    exportStandardThermoModelMineralHKF(m);
    exportStandardThermoModelWaterHKF(m);
    exportStandardThermoModelYAML(m);
}
