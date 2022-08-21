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
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurfaceProps.hpp>

using namespace Reaktoro;

void exportComplexationSurfaceProps(py::module& m)
{
    py::class_<ComplexationSurfaceProps>(m, "ComplexationSurfaceProps")
        .def(py::init<const ComplexationSurface&, const ChemicalSystem&>())
        .def(py::init<const ComplexationSurface&, const ChemicalState&>())
        .def(py::init<const ComplexationSurface&, const ChemicalProps&>())
        .def("update", py::overload_cast<const ChemicalState&>(&ComplexationSurfaceProps::update))
        .def("update", py::overload_cast<const ChemicalProps&>(&ComplexationSurfaceProps::update))
        .def("elementAmounts", &ComplexationSurfaceProps::elementAmounts)
        .def("elementAmount", &ComplexationSurfaceProps::elementAmount)
        .def("speciesAmount", &ComplexationSurfaceProps::speciesAmount)
        .def("speciesAmounts", &ComplexationSurfaceProps::speciesAmounts)
        .def("speciesFractions", &ComplexationSurfaceProps::speciesFractions)
        .def("speciesFraction", &ComplexationSurfaceProps::speciesFraction)
        .def("output", py::overload_cast<std::ostream&>(&ComplexationSurfaceProps::output, py::const_))
        .def("output", py::overload_cast<const String&>(&ComplexationSurfaceProps::output, py::const_))
        .def("__repr__", [](const ComplexationSurfaceProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
