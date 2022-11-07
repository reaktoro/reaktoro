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
#include <Reaktoro/Thermodynamics/Surface/SurfaceSite.hpp>

using namespace Reaktoro;

void exportSurfaceSite(py::module& m)
{
    py::class_<SurfaceSiteState>(m, "SurfaceSiteState")
        .def(py::init<>())
        .def_readwrite("T", &SurfaceSiteState::T)
        .def_readwrite("P", &SurfaceSiteState::P)
        .def_readwrite("x", &SurfaceSiteState::x)
        .def_readwrite("charge", &SurfaceSiteState::charge)
        .def_readwrite("sigma", &SurfaceSiteState::sigma)
        ;
    
    py::class_<SurfaceSite>(m, "SurfaceSite")
        .def(py::init<>())
        .def(py::init<const String&>())
        .def(py::init<const String&, const String&>())
        .def("name", &SurfaceSite::name)
        .def("tag", &SurfaceSite::tag)
        .def("specificSurfaceArea", &SurfaceSite::specificSurfaceArea)
        .def("mass", &SurfaceSite::mass)
        .def("amount", &SurfaceSite::amount)
        .def("species", py::overload_cast<Index>(&SurfaceSite::species, py::const_), return_internal_ref)
        .def("species", py::overload_cast<>(&SurfaceSite::species, py::const_), return_internal_ref)
        .def("speciesIndices", &SurfaceSite::speciesIndices)
        .def("setName", &SurfaceSite::setName)
        .def("setSurfaceName", &SurfaceSite::setSurfaceName)
        .def("setTag", &SurfaceSite::setTag)
        .def("setAmount", &SurfaceSite::setAmount)
        .def("setDensity", &SurfaceSite::setDensity)
        .def("addSorptionSpecies", &SurfaceSite::addSorptionSpecies)
        ;
}
