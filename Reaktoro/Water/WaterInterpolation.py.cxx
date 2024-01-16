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
#include <Reaktoro/Water/WaterInterpolation.hpp>
using namespace Reaktoro;

void exportWaterInterpolation(py::module& m)
{
    m.def("waterDensityWagnerPrussInterp", waterDensityWagnerPrussInterp, "Compute the density of water (in kg/m3) at given a temperature and pressure using quadratic interpolation.");
    m.def("waterThermoPropsWagnerPrussInterp", waterThermoPropsWagnerPrussInterp, "Compute the thermodynamic properties of water at given a temperature and pressure using quadratic interpolation.");
    m.def("waterThermoPropsWagnerPrussInterpData", waterThermoPropsWagnerPrussInterpData, "Return the pre-computed thermodynamic properties of water using Wagner and Pruss (2002) equation of state used for interpolation.");
}
