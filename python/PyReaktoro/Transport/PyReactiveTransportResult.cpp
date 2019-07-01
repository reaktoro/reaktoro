// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Transport/ReactiveTransportResult.hpp>

namespace Reaktoro {

void exportReactiveTransportResult(py::module& m)
{
    py::class_<ReactiveTransportTiming>(m, "ReactiveTransportTiming")
        .def_readwrite("step", &ReactiveTransportTiming::step)
        .def_readwrite("transport", &ReactiveTransportTiming::transport)
        .def_readwrite("equilibrium", &ReactiveTransportTiming::equilibrium)
        .def(py::self += py::self)
        ;

    py::class_<ReactiveTransportResult>(m, "ReactiveTransportResult")
        .def_readwrite("equilibrium_at_cell", &ReactiveTransportResult::equilibrium_at_cell)
        .def_readwrite("smart_equilibrium_at_cell", &ReactiveTransportResult::smart_equilibrium_at_cell)
        .def_readwrite("timing", &ReactiveTransportResult::timing)
        ;
}

} // namespace Reaktoro
