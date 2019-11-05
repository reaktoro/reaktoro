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
#include <Reaktoro/Optimization/OptimumResult.hpp>

namespace Reaktoro {

void exportOptimumResult(py::module& m)
{
    py::class_<OptimumResult>(m, "OptimumResult")
        .def(py::init<>())
        .def_readwrite("succeeded", &OptimumResult::succeeded)
        .def_readwrite("iterations", &OptimumResult::iterations)
        .def_readwrite("num_objective_evals", &OptimumResult::num_objective_evals)
        .def_readwrite("convergence_rate", &OptimumResult::convergence_rate)
        .def_readwrite("error", &OptimumResult::error)
        .def_readwrite("time", &OptimumResult::time)
        .def_readwrite("time_objective_evals", &OptimumResult::time_objective_evals)
        .def_readwrite("time_constraint_evals", &OptimumResult::time_constraint_evals)
        .def_readwrite("time_linear_systems", &OptimumResult::time_linear_systems);
}

} // namespace Reaktoro
