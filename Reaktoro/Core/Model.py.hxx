// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/Model.hpp>
using namespace Reaktoro;

/// Export a Model to Python. This method only exports the methods in class Model. This is useful in
/// case you want to attach custom constructors that work differently from Python.
template<typename Result, typename... Args>
auto exportModelMethodsOnly(py::module& m, const char* modelname)
{
    using ResultRef = Ref<Result>;
    using ModelType = Model<Result(Args...)>;

    return py::class_<ModelType>(m, modelname)
        .def("params", &ModelType::params, return_internal_ref)
        .def("apply", &ModelType::apply)
        .def("__call__", py::overload_cast<const Args&...>(&ModelType::operator(), py::const_))
        .def("__call__", py::overload_cast<ResultRef, const Args&...>(&ModelType::operator(), py::const_))
        ;
}

/// Export a Model to Python, including methods and constructors.
template<typename Result, typename... Args>
auto exportModel(py::module& m, const char* modelname)
{
    using ResultRef = Ref<Result>;
    using ModelType = Model<Result(Args...)>;

    auto cls = exportModelMethodsOnly<Result, Args...>(m, modelname)
        .def(py::init<>())
        .def(py::init<const Fn<Result(Args...)>&>())
        ;

    // At the moment, the other constructor overload is not possible:
    // .def(py::init<const Fn<void(ResultRef, Args...)>&>())

    py::implicitly_convertible<Fn<Result(Args...)>, ModelType>();

    return cls;
}
