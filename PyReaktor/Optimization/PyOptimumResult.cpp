// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "PyOptimumResult.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Optimization/OptimumResult.hpp>

namespace Reaktor {

auto export_OptimumResult() -> void
{
    py::class_<OptimumSolution>("OptimumSolution")
        .def_readwrite("x", &OptimumSolution::x)
        .def_readwrite("y", &OptimumSolution::y)
        .def_readwrite("zl", &OptimumSolution::z)
        .def_readwrite("zu", &OptimumSolution::zu)
        ;

    py::class_<OptimumStatistics>("OptimumStatistics")
        .def_readwrite("converged", &OptimumStatistics::converged)
        .def_readwrite("num_iterations", &OptimumStatistics::num_iterations)
        .def_readwrite("num_objective_evals", &OptimumStatistics::num_objective_evals)
        .def_readwrite("convergence_rate", &OptimumStatistics::convergence_rate)
        .def_readwrite("error", &OptimumStatistics::error)
        .def_readwrite("time", &OptimumStatistics::time)
        .def_readwrite("time_objective_evals", &OptimumStatistics::time_objective_evals)
        .def_readwrite("time_constraint_evals", &OptimumStatistics::time_constraint_evals)
        .def_readwrite("time_linear_system_solutions", &OptimumStatistics::time_linear_system)
        ;

    py::class_<OptimumResult>("OptimumResult")
        .def_readwrite("solution", &OptimumResult::solution)
        .def_readwrite("statistics", &OptimumResult::statistics)
        ;
}

} // namespace Reaktor
