// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Thermodynamics/Water/WaterHelmholtzState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>
using namespace Reaktoro;

void exportWaterThermoStateUtils(py::module& m)
{
    m.def("waterThermoStateHGK", waterThermoStateHGK);
    m.def("waterThermoStateWagnerPruss", waterThermoStateWagnerPruss);
    m.def("waterThermoStateHGKMemoized", waterThermoStateHGKMemoized);
    m.def("waterThermoStateWagnerPrussMemoized", waterThermoStateWagnerPrussMemoized);
    m.def("waterThermoState", waterThermoState);
}
