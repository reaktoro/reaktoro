// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Water/WaterThermoProps.hpp>
using namespace Reaktoro;

void exportWaterThermoProps(py::module& m)
{
    py::class_<WaterThermoProps>(m, "WaterThermoProps")
        .def(py::init<>())
        .def_readwrite("T"  , &WaterThermoProps::T  , "The temperature of water (in K).")
        .def_readwrite("V"  , &WaterThermoProps::V  , "The specific volume of water (in m3/kg).")
        .def_readwrite("S"  , &WaterThermoProps::S  , "The specific entropy of water (in J/(kg*K)).")
        .def_readwrite("A"  , &WaterThermoProps::A  , "The specific Helmholtz free energy of water (in J/kg).")
        .def_readwrite("U"  , &WaterThermoProps::U  , "The specific internal energy of water (in J/kg).")
        .def_readwrite("H"  , &WaterThermoProps::H  , "The specific enthalpy of water (in J/kg).")
        .def_readwrite("G"  , &WaterThermoProps::G  , "The specific Gibbs free energy of water (in J/kg).")
        .def_readwrite("Cv" , &WaterThermoProps::Cv , "The specific isochoric heat capacity of water (in J/(kg*K)).")
        .def_readwrite("Cp" , &WaterThermoProps::Cp , "The specific isobaric heat capacity of water (in J/(kg*K)).")
        .def_readwrite("D"  , &WaterThermoProps::D  , "The specific density of water (in kg/m3).")
        .def_readwrite("DT" , &WaterThermoProps::DT , "The first-order partial derivative of density with respect to temperature (in (kg/m3)/K).")
        .def_readwrite("DP" , &WaterThermoProps::DP , "The first-order partial derivative of density with respect to pressure (in (kg/m3)/Pa).")
        .def_readwrite("DTT", &WaterThermoProps::DTT, "The second-order partial derivative of density with respect to temperature (in (kg/m3)/(K*K)).")
        .def_readwrite("DTP", &WaterThermoProps::DTP, "The second-order partial derivative of density with respect to temperature and pressure (in (kg/m3)/(K*Pa)).")
        .def_readwrite("DPP", &WaterThermoProps::DPP, "The second-order partial derivative of density with respect to pressure (in (kg/m3)/(Pa*Pa)).")
        .def_readwrite("P"  , &WaterThermoProps::P  , "The pressure of water (in Pa).")
        .def_readwrite("PT" , &WaterThermoProps::PT , "The first-order partial derivative of pressure with respect to temperature (in Pa/K).")
        .def_readwrite("PD" , &WaterThermoProps::PD , "The first-order partial derivative of pressure with respect to density (in Pa/(kg/m3)).")
        .def_readwrite("PTT", &WaterThermoProps::PTT, "The second-order partial derivative of pressure with respect to temperature (in Pa/(K*K)).")
        .def_readwrite("PTD", &WaterThermoProps::PTD, "The second-order partial derivative of pressure with respect to temperature and density (in Pa/(K*kg/m3)).")
        .def_readwrite("PDD", &WaterThermoProps::PDD, "The second-order partial derivative of pressure with respect to density (in Pa/((kg/m3)*(kg/m3))).")
        ;
}
