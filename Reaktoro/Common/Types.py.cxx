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
#include <Reaktoro/Common/Types.hpp>
using namespace Reaktoro;

void exportTypes(py::module& m)
{
    const auto __getitem__ = [](const Indices& self, std::size_t i) { return self[i]; };
    const auto __setitem__ = [](Indices& self, std::size_t i, Index val) { self[i] = val; };

    const auto push_back = [](Indices& self, const Index& i) { self.push_back(i); };

    auto createIndices = [](py::list l)
    {
        Indices res(l.size());
        for(auto i = 0; i < res.size(); ++i)
            res[i] = l[i].cast<Index>();
        return res;
    };

    py::class_<Indices>(m, "Indices")
        .def(py::init<>())
        .def(py::init(createIndices))
        .def("empty", &Indices::empty)
        .def("size", &Indices::size)
        .def("pop_back", &Indices::pop_back)
        .def("push_back", push_back)
        .def("append", push_back)
        .def("front", [](const Indices& self) { return self.front(); }, return_internal_ref)
        .def("front", [](Indices& self) { return self.front(); }, return_internal_ref)
        .def("back", [](const Indices& self) { return self.back(); }, return_internal_ref)
        .def("back", [](Indices& self) { return self.back(); }, return_internal_ref)
        .def("__len__", &Indices::size)
        .def("__iter__", [](Indices& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>())
        .def("__getitem__", __getitem__, return_internal_ref)
        .def("__setitem__", __setitem__)
        .def(py::self == py::self)
        .def(py::self != py::self)
        ;

    py::implicitly_convertible<py::list, Indices>();
}
