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
#include <Reaktoro/Core/Surface.hpp>
#include <Reaktoro/Core/SurfaceList.hpp>
using namespace Reaktoro;

void exportSurfaceList(py::module& m)
{
    py::class_<SurfaceList>(m, "SurfaceList")
        .def(py::init<>())
        .def(py::init<const Vec<Surface>&>())
        .def("append", &SurfaceList::append)
        .def("data", &SurfaceList::data)
        .def("empty", &SurfaceList::empty)
        .def("size", &SurfaceList::size)
        .def("find", &SurfaceList::find)
        .def("findWithName", &SurfaceList::findWithName)
        .def("findWithPhases", &SurfaceList::findWithPhases)
        .def("index", &SurfaceList::index)
        .def("indexWithName", &SurfaceList::indexWithName)
        .def("indexWithPhases", &SurfaceList::indexWithPhases)
        .def("get", &SurfaceList::get, return_internal_ref)
        .def("getWithName", &SurfaceList::getWithName, return_internal_ref)
        .def("getWithPhases", &SurfaceList::getWithPhases, return_internal_ref)
        .def("withNames", &SurfaceList::withNames)
        .def("__len__", &SurfaceList::size)
        .def("__getitem__", [](const SurfaceList& self, Index i) { return self[i]; }, return_internal_ref)
        .def("__getitem__", [](SurfaceList& self, Index i) { return self[i]; }, return_internal_ref)
        .def("__iter__", [](const SurfaceList& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()) // keep object alive while iterator exists;
        .def(py::self + py::self);
        ;

    py::implicitly_convertible<Vec<Surface>, SurfaceList>();
}
