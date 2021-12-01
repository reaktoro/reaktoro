// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Transport/TransportResult.hpp>
using namespace Reaktoro;

void exportTransportResult(py::module& m)
{
    py::class_<TransportTiming>(m, "TransportTiming")
        .def_readwrite("step", &TransportTiming::step)
        .def_readwrite("matrix_equation_assembly", &TransportTiming::matrix_equation_assembly)
        .def_readwrite("matrix_equation_solve", &TransportTiming::matrix_equation_solve)
        .def(py::self += py::self)
        ;

    py::class_<TransportResult>(m, "TransportResult")
        .def_readwrite("timing", &TransportResult::timing)
        ;
}
