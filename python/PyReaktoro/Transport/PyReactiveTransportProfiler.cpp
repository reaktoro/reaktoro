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
#include <Reaktoro/Transport/ReactiveTransportProfiler.hpp>
#include <Reaktoro/Transport/ReactiveTransportSolver.hpp>

namespace Reaktoro {

void exportReactiveTransportProfiler(py::module& m)
{
    py::class_<ReactiveTransportProfiler::ComputingCostsPerTimeStep>(m, "ReactiveTransportProfilerComputingCostsPerTimeStep")
        .def_readwrite("t", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::t)
        .def_readwrite("transport", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::transport)
        .def_readwrite("equilibrium", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::equilibrium)
        .def_readwrite("smart_equilibrium", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::smart_equilibrium)
        .def_readwrite("smart_equilibrium_with_ideal_search", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::smart_equilibrium_with_ideal_search)
        .def_readwrite("smart_equilibrium_estimate", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::smart_equilibrium_estimate)
        .def_readwrite("smart_equilibrium_nearest_neighbor_search", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::smart_equilibrium_nearest_neighbor_search)
        .def_readwrite("smart_equilibrium_gibbs_energy_minimization", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::smart_equilibrium_gibbs_energy_minimization)
        .def_readwrite("smart_equilibrium_storage", &ReactiveTransportProfiler::ComputingCostsPerTimeStep::smart_equilibrium_storage)
        ;

    py::class_<ReactiveTransportProfiler::SmartEquilibriumProfiling>(m, "ReactiveTransportProfilerSmartEquilibriumProfiling")
        .def_readwrite("timing", &ReactiveTransportProfiler::SmartEquilibriumProfiling::timing)
        .def_readwrite("num_equilibrium_calculations", &ReactiveTransportProfiler::SmartEquilibriumProfiling::num_equilibrium_calculations)
        .def_readwrite("num_smart_equilibrium_accepted_estimates", &ReactiveTransportProfiler::SmartEquilibriumProfiling::num_smart_equilibrium_accepted_estimates)
        .def_readwrite("num_smart_equilibrium_required_learnings", &ReactiveTransportProfiler::SmartEquilibriumProfiling::num_smart_equilibrium_required_learnings)
        .def_readwrite("smart_equilibrium_estimate_acceptance_rate", &ReactiveTransportProfiler::SmartEquilibriumProfiling::smart_equilibrium_estimate_acceptance_rate)
        .def_readwrite("cells_where_learning_was_required_at_step", &ReactiveTransportProfiler::SmartEquilibriumProfiling::cells_where_learning_was_required_at_step)
        ;

    py::class_<ReactiveTransportProfiler>(m, "ReactiveTransportProfiler")
        .def(py::init<const ReactiveTransportSolver&>())
        .def("update", &ReactiveTransportProfiler::update)
        .def("computingCostsPerTimeStep", &ReactiveTransportProfiler::computingCostsPerTimeStep)
        .def("smartEquilibriumProfiling", &ReactiveTransportProfiler::smartEquilibriumProfiling)
        .def("results", &ReactiveTransportProfiler::results)
        ;
}

} // namespace Reaktoro
