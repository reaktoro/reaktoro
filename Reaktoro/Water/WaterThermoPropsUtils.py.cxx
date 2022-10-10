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

// Reaktoro includes
#include <Reaktoro/Water/WaterHelmholtzProps.hpp>
#include <Reaktoro/Water/WaterThermoProps.hpp>
#include <Reaktoro/Water/WaterThermoPropsUtils.hpp>
using namespace Reaktoro;

void exportWaterThermoPropsUtils(py::module& m)
{
    m.def("waterThermoPropsHGK", waterThermoPropsHGK, "Calculate the thermodynamic properties of water using the Haar-Gallagher-Kell (1984) equation of state.");
    m.def("waterThermoPropsWagnerPruss", waterThermoPropsWagnerPruss, "Calculate the thermodynamic properties of water using the Haar-Gallagher-Kell (1984) equation of state.");
    m.def("waterThermoPropsHGKMemoized", waterThermoPropsHGKMemoized, "Calculate the thermodynamic properties of water using the Wagner and Pruss (1995) equation of state.");
    m.def("waterThermoPropsWagnerPrussMemoized", waterThermoPropsWagnerPrussMemoized, "Calculate the thermodynamic properties of water using the Wagner and Pruss (1995) equation of state.");
    m.def("waterThermoPropsWagnerPrussInterpMemoized", waterThermoPropsWagnerPrussInterpMemoized, "Calculate the thermodynamic properties of water using interpolation of pre-computed properties using the Wagner and Pruss (1995) equation of state.");
    m.def("waterThermoProps", waterThermoProps, "Calculate the thermodynamic properties of water.");
}
