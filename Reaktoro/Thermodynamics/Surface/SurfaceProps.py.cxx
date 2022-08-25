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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Thermodynamics/Surface/Surface.hpp>
#include <Reaktoro/Thermodynamics/Surface/SurfaceProps.hpp>

using namespace Reaktoro;

void exportSurfaceProps(py::module& m)
{
    py::class_<SurfaceProps>(m, "SurfaceProps")
        .def(py::init<const Surface&, const ChemicalSystem&>())
        .def(py::init<const Surface&, const ChemicalState&>())
        .def("update", py::overload_cast<const ChemicalState&>(&SurfaceProps::update))
        .def("elementAmounts", &SurfaceProps::elementAmounts)
        .def("elementAmount", &SurfaceProps::elementAmount)
        .def("speciesAmount", &SurfaceProps::speciesAmount)
        .def("speciesAmounts", &SurfaceProps::speciesAmounts)
        .def("speciesFractions", &SurfaceProps::speciesFractions)
        .def("speciesFraction", &SurfaceProps::speciesFraction)
        .def("output", py::overload_cast<std::ostream&>(&SurfaceProps::output, py::const_))
        .def("output", py::overload_cast<const String&>(&SurfaceProps::output, py::const_))
        .def("__repr__", [](const SurfaceProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
