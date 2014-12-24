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

#include "PyEquilibriumResult.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

namespace Reaktor {

auto export_EquilibriumResult() -> void
{
	py::class_<EquilibriumSolution>("EquilibriumSolution")
        .def_readwrite("n", &EquilibriumSolution::n)
        .def_readwrite("y", &EquilibriumSolution::y)
        .def_readwrite("z", &EquilibriumSolution::z)
        .def_readwrite("dndt", &EquilibriumSolution::dndt)
        .def_readwrite("dndp", &EquilibriumSolution::dndp)
        .def_readwrite("dndb", &EquilibriumSolution::dndb)
        .def_readwrite("indices_stable_species", &EquilibriumSolution::indices_stable_species)
		;

    py::class_<EquilibriumStatistics, py::bases<OptimumStatistics>>("EquilibriumStatistics")
        .def(py::init<>())
        .def(py::init<const OptimumStatistics&>())
        ;

    py::class_<EquilibriumResult>("EquilibriumResult")
        .def_readwrite("solution", &EquilibriumResult::solution)
        .def_readwrite("statistics", &EquilibriumResult::statistics)
        ;
}

} // namespace Reaktor
