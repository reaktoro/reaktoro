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

// Models/ActivityModels/Support
void exportAqueousMixture(py::module& m);
void exportCubicEOS(py::module& m);
void exportIonExchangeSurface(py::module& m);

// Models/ActivityModels
void exportActivityModelCubicEOS(py::module& m);
void exportActivityModelDebyeHuckel(py::module& m);
void exportActivityModelDrummond(py::module& m);
void exportActivityModelDuanSun(py::module& m);
void exportActivityModelHKF(py::module& m);
void exportActivityModelIdealAqueous(py::module& m);
void exportActivityModelIdealGas(py::module& m);
void exportActivityModelIdealIonExchange(py::module& m);
void exportActivityModelIdealSolution(py::module& m);
void exportActivityModelIonExchange(py::module& m);
void exportActivityModelPitzerHMW(py::module& m);
void exportActivityModelRedlichKister(py::module& m);
void exportActivityModelRumpf(py::module& m);
void exportActivityModelSetschenow(py::module& m);
void exportActivityModelSpycherPruessEnnis(py::module& m);
void exportActivityModelSpycherReed(py::module& m);
void exportActivityModelVanLaar(py::module& m);

// Models/ReactionThermoModels/Support
void exportReactionThermoModelConstLgK(py::module& m);
void exportReactionThermoModelGemsLgK(py::module& m);
void exportReactionThermoModelPhreeqcLgK(py::module& m);
void exportReactionThermoModelPressureCorrection(py::module& m);
void exportReactionThermoModelVantHoff(py::module& m);
void exportReactionThermoModelYAML(py::module& m);

// Models/StandardThermoModels/Support
void exportStandardThermoModelConstant(py::module& m);
void exportStandardThermoModelHKF(py::module& m);
void exportStandardThermoModelHollandPowell(py::module& m);
void exportStandardThermoModelInterpolation(py::module& m);
void exportStandardThermoModelMaierKelley(py::module& m);
void exportStandardThermoModelMineralHKF(py::module& m);
void exportStandardThermoModelNasa(py::module& m);
void exportStandardThermoModelWaterHKF(py::module& m);
void exportStandardThermoModelYAML(py::module& m);

void exportModels(py::module& m)
{
    // Models/ActivityModels/Support
    exportAqueousMixture(m);
    exportCubicEOS(m);
    exportIonExchangeSurface(m);

    // Models/ActivityModels
    exportActivityModelCubicEOS(m);
    exportActivityModelDebyeHuckel(m);
    exportActivityModelDrummond(m);
    exportActivityModelDuanSun(m);
    exportActivityModelHKF(m);
    exportActivityModelIdealAqueous(m);
    exportActivityModelIdealGas(m);
    exportActivityModelIdealIonExchange(m);
    exportActivityModelIdealSolution(m);
    exportActivityModelIonExchange(m);
    exportActivityModelPitzerHMW(m);
    exportActivityModelRedlichKister(m);
    exportActivityModelRumpf(m);
    exportActivityModelSetschenow(m);
    exportActivityModelSpycherPruessEnnis(m);
    exportActivityModelSpycherReed(m);
    exportActivityModelVanLaar(m);

    // Models/ReactionThermoModels/Support
    exportReactionThermoModelConstLgK(m);
    exportReactionThermoModelGemsLgK(m);
    exportReactionThermoModelPhreeqcLgK(m);
    exportReactionThermoModelPressureCorrection(m);
    exportReactionThermoModelVantHoff(m);
    exportReactionThermoModelYAML(m);

    // Models/StandardThermoModels/Support
    exportStandardThermoModelConstant(m);
    exportStandardThermoModelHKF(m);
    exportStandardThermoModelHollandPowell(m);
    exportStandardThermoModelInterpolation(m);
    exportStandardThermoModelMaierKelley(m);
    exportStandardThermoModelMineralHKF(m);
    exportStandardThermoModelNasa(m);
    exportStandardThermoModelWaterHKF(m);
    exportStandardThermoModelYAML(m);
}
