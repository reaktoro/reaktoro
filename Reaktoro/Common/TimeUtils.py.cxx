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
#include <pybind11/chrono.h>

// Reaktoro includes
#include <Reaktoro/Common/TimeUtils.hpp>
using namespace Reaktoro;

void exportTimeUtils(py::module& m)
{
    m.def("time", Reaktoro::time);
    m.def("elapsed", py::overload_cast<Time const&, Time const&>(elapsed));
    m.def("elapsed", py::overload_cast<Time const&>(elapsed));

    py::class_<Stopwatch>(m, "Stopwatch")
        .def(py::init<>())
        .def("start", &Stopwatch::start, "Start measuring time.")
        .def("pause", &Stopwatch::pause, "Pause measuring time.")
        .def("reset", &Stopwatch::reset, "Reset the stopwatch.")
        .def("time", &Stopwatch::time, "Get the accumulated elapsed time (in seconds) between calls to methods start and pause.")
        ;
}
