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

void exportStandardThermoModelConstant(py::module& m);
void exportStandardThermoModelHKF(py::module& m);
void exportStandardThermoModelHollandPowell(py::module& m);
void exportStandardThermoModelInterpolation(py::module& m);
void exportStandardThermoModelMaierKelley(py::module& m);
void exportStandardThermoModelMineralHKF(py::module& m);
void exportStandardThermoModelNasa(py::module& m);
void exportStandardThermoModelWaterHKF(py::module& m);
void exportStandardThermoModelFromData(py::module& m);

void exportReactionStandardThermoModelConstLgK(py::module& m);
void exportReactionStandardThermoModelGemsLgK(py::module& m);
void exportReactionStandardThermoModelPhreeqcLgK(py::module& m);
void exportReactionStandardThermoModelPressureCorrection(py::module& m);
void exportReactionStandardThermoModelVantHoff(py::module& m);
void exportReactionStandardThermoModelFromData(py::module& m);

void exportStandardVolumeModelConstant(py::module& m);

void exportStandardThermoModels(py::module& m)
{
    exportStandardThermoModelConstant(m);
    exportStandardThermoModelHKF(m);
    exportStandardThermoModelHollandPowell(m);
    exportStandardThermoModelInterpolation(m);
    exportStandardThermoModelMaierKelley(m);
    exportStandardThermoModelMineralHKF(m);
    exportStandardThermoModelNasa(m);
    exportStandardThermoModelWaterHKF(m);
    exportStandardThermoModelFromData(m);

    exportReactionStandardThermoModelConstLgK(m);
    exportReactionStandardThermoModelGemsLgK(m);
    exportReactionStandardThermoModelPhreeqcLgK(m);
    exportReactionStandardThermoModelPressureCorrection(m);
    exportReactionStandardThermoModelVantHoff(m);
    exportReactionStandardThermoModelFromData(m);

    exportStandardVolumeModelConstant(m);
}
