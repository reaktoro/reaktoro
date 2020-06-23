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
#include <Reaktoro/Transport/ReactiveTransportAnalysis.hpp>

namespace Reaktoro {

void exportReactiveTransportAnalysis(py::module& m)
{
    py::class_<ReactiveTransportAnalysis>(m, "ReactiveTransportAnalysis")
            .def(py::init<>())
            .def_readwrite("transport", &ReactiveTransportAnalysis::transport)
            .def_readwrite("equilibrium", &ReactiveTransportAnalysis::equilibrium)
            .def_readwrite("smart_equilibrium", &ReactiveTransportAnalysis::smart_equilibrium)
            .def_readwrite("computing_costs_per_time_step", &ReactiveTransportAnalysis::computing_costs_per_time_step)
            ;

    py::class_<ReactiveTransportAnalysis::TransportAnalysis>(m, "ReactiveTransportAnalysis_TransportAnalysis")
        .def_readwrite("timing", &ReactiveTransportAnalysis::TransportAnalysis::timing)
        ;

    py::class_<ReactiveTransportAnalysis::SmartEquilibriumAnalysis>(m, "ReactiveTransportAnalysis_SmartEquilibriumAnalysis")
        .def_readwrite("timing", &ReactiveTransportAnalysis::SmartEquilibriumAnalysis::timing)
        .def_readwrite("num_equilibrium_calculations", &ReactiveTransportAnalysis::SmartEquilibriumAnalysis::num_equilibrium_calculations)
        .def_readwrite("num_smart_equilibrium_accepted_estimates", &ReactiveTransportAnalysis::SmartEquilibriumAnalysis::num_smart_equilibrium_accepted_estimates)
        .def_readwrite("num_smart_equilibrium_required_learnings", &ReactiveTransportAnalysis::SmartEquilibriumAnalysis::num_smart_equilibrium_required_learnings)
        .def_readwrite("smart_equilibrium_estimate_acceptance_rate", &ReactiveTransportAnalysis::SmartEquilibriumAnalysis::smart_equilibrium_estimate_acceptance_rate)
        .def_readwrite("cells_where_learning_was_required_at_step", &ReactiveTransportAnalysis::SmartEquilibriumAnalysis::cells_where_learning_was_required_at_step)
        ;

    py::class_<ReactiveTransportAnalysis::ComputingCostsPerTimeStep>(m, "ReactiveTransportAnalysis_ComputingCostsPerTimeStep")
        .def_readwrite("transport", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::transport)
        .def_readwrite("equilibrium", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::equilibrium)
        .def_readwrite("smart_equilibrium", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::smart_equilibrium)
        .def_readwrite("smart_equilibrium_with_ideal_search", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::smart_equilibrium_with_ideal_search)
        .def_readwrite("smart_equilibrium_estimate", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::smart_equilibrium_estimate)
        .def_readwrite("smart_equilibrium_search", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::smart_equilibrium_search)
        .def_readwrite("smart_equilibrium_gibbs_energy_minimization", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::smart_equilibrium_gibbs_energy_minimization)
        .def_readwrite("smart_equilibrium_storage", &ReactiveTransportAnalysis::ComputingCostsPerTimeStep::smart_equilibrium_storage)
        ;
}

} // namespace Reaktoro
