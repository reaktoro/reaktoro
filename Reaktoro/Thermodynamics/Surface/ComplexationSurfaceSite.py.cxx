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
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurfaceSite.hpp>

using namespace Reaktoro;

void exportComplexationSurfaceSite(py::module& m)
{
    py::class_<ComplexationSurfaceSiteState>(m, "ComplexationSurfaceSiteState")
        .def(py::init<>())
        .def_readwrite("T", &ComplexationSurfaceSiteState::T)
        .def_readwrite("P", &ComplexationSurfaceSiteState::P)
        .def_readwrite("x", &ComplexationSurfaceSiteState::x)
        .def_readwrite("charge", &ComplexationSurfaceSiteState::charge)
        .def_readwrite("sigma", &ComplexationSurfaceSiteState::sigma)
        ;
    
    py::class_<ComplexationSurfaceSite>(m, "ComplexationSurfaceSite")
        .def(py::init<>())
        .def(py::init<const String&>())
        .def(py::init<const String&, const String&>())
        .def("name", &ComplexationSurfaceSite::name)
        .def("tag", &ComplexationSurfaceSite::tag)
        .def("specificSurfaceArea", &ComplexationSurfaceSite::specificSurfaceArea)
        .def("mass", &ComplexationSurfaceSite::mass)
        .def("amount", &ComplexationSurfaceSite::amount)
        .def("species", py::overload_cast<Index>(&ComplexationSurfaceSite::species, py::const_), return_internal_ref)
        .def("species", py::overload_cast<>(&ComplexationSurfaceSite::species, py::const_), return_internal_ref)
        .def("speciesIndices", &ComplexationSurfaceSite::speciesIndices)
        .def("setName", &ComplexationSurfaceSite::setName)
        .def("setSurfaceName", &ComplexationSurfaceSite::setSurfaceName)
        .def("setTag", &ComplexationSurfaceSite::setTag)
        .def("setAmount", &ComplexationSurfaceSite::setAmount)
        .def("setDensity", &ComplexationSurfaceSite::setDensity)
        .def("addSorptionSpecies", &ComplexationSurfaceSite::addSorptionSpecies)
        ;
}
