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

#include "PyEquilibriumSolver.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktor/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>
#include <Reaktor/Equilibrium/EquilibriumSolver.hpp>

namespace Reaktor {

auto export_EquilibriumSolver() -> void
{
    void (*initialize1)(const EquilibriumProblem&, EquilibriumResult&) = initialize;
    void (*initialize2)(const EquilibriumProblem&, EquilibriumResult&, const EquilibriumOptions&) = initialize;
    void (*solve1)(const EquilibriumProblem&, EquilibriumResult&) = solve;
    void (*solve2)(const EquilibriumProblem&, EquilibriumResult&, const EquilibriumOptions&) = solve;
    
    py::def("initialize", initialize1);
    py::def("initialize", initialize2);
    py::def("solve", solve1);
    py::def("solve", solve2);
}

} // namespace Reaktor
