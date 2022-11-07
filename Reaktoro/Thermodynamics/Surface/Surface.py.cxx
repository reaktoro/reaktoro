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
#include <Reaktoro/Thermodynamics/Surface/Surface.hpp>
#include <Reaktoro/Thermodynamics/Surface/SurfaceSite.hpp>

using namespace Reaktoro;

void exportSurface(py::module& m)
{
    py::class_<SurfaceState>(m, "SurfaceState")
        .def(py::init<>())
        .def("updatePotential", &SurfaceState::updatePotential)
        .def_readwrite("T", &SurfaceState::T)
        .def_readwrite("P", &SurfaceState::P)
        .def_readwrite("x", &SurfaceState::x)
        .def_readwrite("charge", &SurfaceState::charge)
        .def_readwrite("sigma", &SurfaceState::sigma)
        .def_readwrite("As", &SurfaceState::As)
        .def_readwrite("mass", &SurfaceState::mass)
        .def_readwrite("psi", &SurfaceState::psi)
        ;
    py::class_<Surface>(m, "Surface")
        .def(py::init<const String&>())
        .def(py::init<const SpeciesList&>())
        .def("clone", &Surface::clone, return_internal_ref)
        .def("name", &Surface::name)
        .def("potential", py::overload_cast<>(&Surface::potential, py::const_), return_internal_ref)
        .def("potential", py::overload_cast<real, real, real>(&Surface::potential, py::const_), return_internal_ref)
        .def("species", py::overload_cast<Index>(&Surface::species, py::const_), return_internal_ref)
        .def("species", py::overload_cast<>(&Surface::species, py::const_), return_internal_ref)
        .def("charges", &Surface::charges)
        .def("moleFractions", &Surface::moleFractions)
        .def("specificSurfaceArea", &Surface::specificSurfaceArea)
        .def("mass", &Surface::mass)
        .def("sites", &Surface::sites)
        .def("addSurfaceSpecies", &Surface::addSurfaceSpecies)
        .def("setName", &Surface::setName)
        .def("setSpecificSurfaceArea", &Surface::setSpecificSurfaceArea)
        .def("setMass", &Surface::setMass)
        .def("addSite", py::overload_cast<const String&, const String&>(&Surface::addSite), return_internal_ref)
        .def("addSite", py::overload_cast<const SurfaceSite&>(&Surface::addSite), return_internal_ref)
        .def("__repr__", [](const Surface& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
