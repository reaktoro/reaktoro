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
#include <Reaktoro/Transport/ReactiveTransportAnalysis.hpp>
using namespace Reaktoro;

void exportReactiveTransportAnalysis(py::module& m)
{
    py::class_<ReactiveTransportAnalysis>(m, "ReactiveTransportAnalysis")
        .def(py::init<>())
        .def_readwrite("transport", &ReactiveTransportAnalysis::transport)
        .def_readwrite("equilibrium", &ReactiveTransportAnalysis::equilibrium)
        .def_readwrite("computing_costs_per_time_step", &ReactiveTransportAnalysis::computing_costs_per_time_step)
        ;

    py::class_<ReactiveTransportAnalysis::TransportAnalysis>(m, "ReactiveTransportAnalysis_TransportAnalysis")
        .def_readwrite("timing", &ReactiveTransportAnalysis::TransportAnalysis::timing)
        ;

    py::class_<ReactiveTransportAnalysis::ComputingCostsPerTimeStep>(m, "ReactiveTransportAnalysis_ComputingCostsPerTimeStep")
        .def_readwrite("transport", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::transport)
        .def_readwrite("equilibrium", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::equilibrium)
        ;
}