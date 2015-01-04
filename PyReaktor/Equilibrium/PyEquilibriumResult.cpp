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
#include <Reaktor/Equilibrium/EquilibriumResult.hpp>

namespace Reaktor {

auto export_EquilibriumResult() -> void
{
	py::class_<EquilibriumState>("EquilibriumState")
        .def_readwrite("n", &EquilibriumState::n)
        .def_readwrite("y", &EquilibriumState::ye)
        .def_readwrite("z", &EquilibriumState::ze)
        .def_readwrite("dndt", &EquilibriumState::dndt)
        .def_readwrite("dndp", &EquilibriumState::dndp)
        .def_readwrite("dndb", &EquilibriumState::dndb)
        .def_readwrite("indices_stable_species", &EquilibriumState::indices_stable_species)
		;

    py::class_<EquilibriumResult, py::bases<OptimumResult>>("EquilibriumStatistics")
        .def(py::init<>())
        .def(py::init<const OptimumResult&>())
        ;

    py::class_<EquilibriumState>("EquilibriumResult")
        .def_readwrite("solution", &EquilibriumState::solution)
        .def_readwrite("statistics", &EquilibriumState::statistics)
        ;
}

} // namespace Reaktor
