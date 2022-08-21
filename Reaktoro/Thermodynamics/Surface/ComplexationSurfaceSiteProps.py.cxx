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
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurfaceSite.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurfaceSiteProps.hpp>

using namespace Reaktoro;

void exportComplexationSurfaceSiteProps(py::module& m)
{
    py::class_<ComplexationSurfaceSiteProps>(m, "ComplexationSurfaceSiteProps")
        .def(py::init<const ComplexationSurfaceSite&, const ChemicalSystem&>())
        .def(py::init<const ComplexationSurfaceSite&, const ChemicalState&>())
        .def(py::init<const ComplexationSurfaceSite&, const ChemicalProps&>())
        .def("update", py::overload_cast<const ChemicalState&>(&ComplexationSurfaceSiteProps::update))
        .def("update", py::overload_cast<const ChemicalProps&>(&ComplexationSurfaceSiteProps::update))
        .def("elementAmounts", &ComplexationSurfaceSiteProps::elementAmounts)
        .def("elementAmount", &ComplexationSurfaceSiteProps::elementAmount)
        .def("speciesAmount", &ComplexationSurfaceSiteProps::speciesAmount)
        .def("speciesAmounts", &ComplexationSurfaceSiteProps::speciesAmounts)
        .def("speciesFractions", &ComplexationSurfaceSiteProps::speciesFractions)
        .def("speciesFraction", &ComplexationSurfaceSiteProps::speciesFraction)
        .def("charge", &ComplexationSurfaceSiteProps::charge)
        .def("sigma", &ComplexationSurfaceSiteProps::sigma)
        .def("output", py::overload_cast<std::ostream&>(&ComplexationSurfaceSiteProps::output, py::const_))
        .def("output", py::overload_cast<const String&>(&ComplexationSurfaceSiteProps::output, py::const_))
        .def("__repr__", [](const ComplexationSurfaceSiteProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
