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

// Thermodynamics/Aqueous
void exportActivityModelDuanSun(py::module& m);
void exportActivityModelPitzerHMW(py::module& m);
void exportAqueousProps(py::module& m);
void exportActivityModelSetschenow(py::module& m);
void exportActivityModelHKF(py::module& m);
void exportActivityModelDrummond(py::module& m);
void exportActivityModelDavies(py::module& m);
void exportActivityModelDebyeHuckel(py::module& m);
void exportActivityModelRumpf(py::module& m);
void exportAqueousMixture(py::module& m);

// Thermodynamics/Fluids
void exportCubicEOS(py::module& m);
void exportActivityModelCubicEOS(py::module& m);
void exportActivityModelSpycherPruessEnnis(py::module& m);
void exportActivityModelSpycherReed(py::module& m);

// Thermodynamics/Ideal
void exportActivityModelIdealAqueous(py::module& m);
void exportActivityModelIdealGas(py::module& m);
void exportActivityModelIdealSolution(py::module& m);
void exportActivityModelIdealIonExchange(py::module& m);

// Thermodynamics/Solids
void exportActivityModelRedlichKister(py::module& m);
void exportActivityModelVanLaar(py::module& m);

// Thermodynamics/Water
void exportWaterConstants(py::module& m);
void exportWaterElectroProps(py::module& m);
void exportWaterElectroPropsJohnsonNorton(py::module& m);
void exportWaterHelmholtzProps(py::module& m);
void exportWaterHelmholtzPropsHGK(py::module& m);
void exportWaterHelmholtzPropsWagnerPruss(py::module& m);
void exportWaterThermoProps(py::module& m);
void exportWaterThermoPropsUtils(py::module& m);
void exportWaterUtils(py::module& m);

// Thermodynamics/Surface
void exportActivityModelIonExchange(py::module& m);
void exportIonExchangeSurface(py::module& m);
void exportIonExchangeProps(py::module& m);


void exportThermodynamics(py::module& m)
{
    // Thermodynamics/Aqueous
    exportActivityModelDuanSun(m);
    exportActivityModelPitzerHMW(m);
    exportAqueousProps(m);
    exportActivityModelSetschenow(m);
    exportActivityModelHKF(m);
    exportActivityModelDrummond(m);
    exportActivityModelDavies(m);
    exportActivityModelDebyeHuckel(m);
    exportActivityModelRumpf(m);
    exportAqueousMixture(m);

    // Thermodynamics/Fluids
    exportCubicEOS(m);
    exportActivityModelCubicEOS(m);
    exportActivityModelSpycherPruessEnnis(m);
    exportActivityModelSpycherReed(m);

    // Thermodynamics/Ideal
    exportActivityModelIdealAqueous(m);
    exportActivityModelIdealGas(m);
    exportActivityModelIdealSolution(m);
    exportActivityModelIdealIonExchange(m);

    // Thermodynamics/Solids
    exportActivityModelRedlichKister(m);
    exportActivityModelVanLaar(m);

    // Thermodynamics/Water
    exportWaterConstants(m);
    exportWaterElectroProps(m);
    exportWaterElectroPropsJohnsonNorton(m);
    exportWaterHelmholtzProps(m);
    exportWaterHelmholtzPropsHGK(m);
    exportWaterHelmholtzPropsWagnerPruss(m);
    exportWaterThermoProps(m);
    exportWaterThermoPropsUtils(m);
    exportWaterUtils(m);

    // Thermodynamics/Surface
    exportActivityModelIonExchange(m);
    exportIonExchangeSurface(m);
    exportIonExchangeProps(m);
}
