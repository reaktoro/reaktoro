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
#include <Reaktoro/Core/Material.hpp>
using namespace Reaktoro;

void exportChemicalFormula(py::module& m)
{
    auto equivalent1 = [](const Material& self, const Material& other)
    {
        return self.equivalent(other);
    };

    auto equivalent2 = [](const Material& l, const Material& r)
    {
        return Material::equivalent(l, r);
    };

    py::class_<Material>(m, "Material")
        .def(py::init<>())
        .def(py::init<String>())
        .def(py::init<String, Pairs<String, double>, double>(),
            py::arg("formula"),
            py::arg("symbols"),
            py::arg("charge"))
        .def("str", &Material::str)
        .def("elements", &Material::elements)
        .def("symbols", &Material::symbols)
        .def("coefficients", &Material::coefficients)
        .def("coefficient", &Material::coefficient)
        .def("charge", &Material::charge)
        .def("equivalent", equivalent1)  // pybind11 does not support overloading both static and instance methods
        // .def_static("equivalent", equivalent2)
        .def(py::self == py::self)
        .def("__repr__", [](const Material& self) { return self.str(); })
        ;

    py::implicitly_convertible<String, Material>();
}
