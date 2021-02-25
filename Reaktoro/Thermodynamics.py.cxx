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
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Thermodynamics/Aqueous
void exportActivityModelDuanSun(py::module& m);
void exportActivityModelPitzerHMW(py::module& m);
void exportAqueousProps(py::module& m);
void exportActivityModelSetschenow(py::module& m);
void exportActivityModelHKF(py::module& m);
void exportActivityModelDrummond(py::module& m);
void exportActivityModelDebyeHuckel(py::module& m);
void exportActivityModelRumpf(py::module& m);
void exportAqueousMixture(py::module& m);

// Thermodynamics/Fluids
void exportActivityModelCubicEOS(py::module& m);
void exportActivityModelSpycherPruessEnnis(py::module& m);
void exportActivityModelSpycherReed(py::module& m);

// Thermodynamics/Ideal
void exportActivityModelIdealAqueous(py::module& m);
void exportActivityModelIdealGas(py::module& m);
void exportActivityModelIdealSolution(py::module& m);

// Thermodynamics/Reactions
void exportReactionThermoModelAnalyticalGEMS(py::module& m);
void exportReactionThermoModelAnalyticalPHREEQC(py::module& m);
void exportReactionThermoModelConstLgK(py::module& m);
void exportReactionThermoModelPressureCorrection(py::module& m);
void exportReactionThermoModelVantHoff(py::module& m);

// Thermodynamics/Solids
void exportActivityModelRedlichKister(py::module& m);
void exportActivityModelVanLaar(py::module& m);

// Thermodynamics/Water
void exportWaterConstants(py::module& m);
void exportWaterElectroState(py::module& m);
void exportWaterElectroStateJohnsonNorton(py::module& m);
void exportWaterHelmholtzState(py::module& m);
void exportWaterHelmholtzStateHGK(py::module& m);
void exportWaterHelmholtzStateWagnerPruss(py::module& m);
void exportWaterThermoState(py::module& m);
void exportWaterThermoStateUtils(py::module& m);
void exportWaterUtils(py::module& m);

void exportThermodynamics(py::module& m)
{
    // Thermodynamics/Aqueous
    exportActivityModelDuanSun(m);
    exportActivityModelPitzerHMW(m);
    exportAqueousProps(m);
    exportActivityModelSetschenow(m);
    exportActivityModelHKF(m);
    exportActivityModelDrummond(m);
    exportActivityModelDebyeHuckel(m);
    exportActivityModelRumpf(m);
    exportAqueousMixture(m);

    // Thermodynamics/Fluids
    exportActivityModelCubicEOS(m);
    exportActivityModelSpycherPruessEnnis(m);
    exportActivityModelSpycherReed(m);

    // Thermodynamics/Ideal
    exportActivityModelIdealAqueous(m);
    exportActivityModelIdealGas(m);
    exportActivityModelIdealSolution(m);

    // Thermodynamics/Reactions
    exportReactionThermoModelAnalyticalGEMS(m);
    exportReactionThermoModelAnalyticalPHREEQC(m);
    exportReactionThermoModelConstLgK(m);
    exportReactionThermoModelPressureCorrection(m);
    exportReactionThermoModelVantHoff(m);

    // Thermodynamics/Solids
    exportActivityModelRedlichKister(m);
    exportActivityModelVanLaar(m);

    // Thermodynamics/Water
    exportWaterConstants(m);
    exportWaterElectroState(m);
    exportWaterElectroStateJohnsonNorton(m);
    exportWaterHelmholtzState(m);
    exportWaterHelmholtzStateHGK(m);
    exportWaterHelmholtzStateWagnerPruss(m);
    exportWaterThermoState(m);
    exportWaterThermoStateUtils(m);
    exportWaterUtils(m);
}
