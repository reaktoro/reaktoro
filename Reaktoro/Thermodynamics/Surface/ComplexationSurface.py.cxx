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
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurfaceSite.hpp>

using namespace Reaktoro;

void exportComplexationSurface(py::module& m)
{
    py::class_<ComplexationSurfaceState>(m, "ComplexationSurfaceState")
        .def(py::init<>())
        .def("updatePotential", &ComplexationSurfaceState::updatePotential)
        .def_readwrite("T", &ComplexationSurfaceState::T)
        .def_readwrite("P", &ComplexationSurfaceState::P)
        .def_readwrite("x", &ComplexationSurfaceState::x)
        .def_readwrite("z", &ComplexationSurfaceState::z)
        .def_readwrite("charge", &ComplexationSurfaceState::charge)
        .def_readwrite("sigma", &ComplexationSurfaceState::sigma)
        .def_readwrite("As", &ComplexationSurfaceState::As)
        .def_readwrite("mass", &ComplexationSurfaceState::mass)
        .def_readwrite("psi", &ComplexationSurfaceState::psi)
        ;
    py::class_<ComplexationSurface>(m, "ComplexationSurface")
        .def(py::init<const String&>())
        .def(py::init<const SpeciesList&>())
        .def("clone", &ComplexationSurface::clone, py::return_value_policy::reference_internal)
        .def("name", &ComplexationSurface::name)
        .def("potential", &ComplexationSurface::potential)
        .def("species", py::overload_cast<Index>(&ComplexationSurface::species, py::const_), py::return_value_policy::reference_internal)
        .def("species", py::overload_cast<>(&ComplexationSurface::species, py::const_), py::return_value_policy::reference_internal)
        .def("charges", &ComplexationSurface::charges)
        .def("equivalentsNumbers", &ComplexationSurface::equivalentsNumbers)
        .def("moleFractions", &ComplexationSurface::moleFractions)
        .def("specificSurfaceArea", &ComplexationSurface::specificSurfaceArea)
        .def("mass", &ComplexationSurface::mass)
        .def("sites", &ComplexationSurface::sites)
        .def("addSurfaceSpecies", &ComplexationSurface::addSurfaceSpecies)
        .def("setName", &ComplexationSurface::setName)
        .def("setMineral", &ComplexationSurface::setMineral)
        .def("setSpecificSurfaceArea", &ComplexationSurface::setSpecificSurfaceArea)
        .def("setMass", &ComplexationSurface::setMass)
        .def("addSite", py::overload_cast<const String&, const String&>(&ComplexationSurface::addSite), py::return_value_policy::reference_internal)
        .def("addSite", py::overload_cast<const ComplexationSurfaceSite&>(&ComplexationSurface::addSite), py::return_value_policy::reference_internal)
        .def("__repr__", [](const ComplexationSurface& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
