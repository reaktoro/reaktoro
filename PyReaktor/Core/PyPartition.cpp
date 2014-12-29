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

#include "PyPartition.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>

namespace Reaktor {

auto export_Partition() -> void
{
	py::class_<Partition>("Partition")
		.def(py::init<>())
		.def(py::init<const Indices&, const Indices&, const Indices&>())
		.def(py::init<const Partition&>())
        .def("equilibriumSpeciesIndices", &Partition::equilibriumSpeciesIndices, py::return_value_policy<py::copy_const_reference>())
        .def("kineticSpeciesIndices", &Partition::kineticSpeciesIndices, py::return_value_policy<py::copy_const_reference>())
        .def("inertSpeciesIndices", &Partition::inertSpeciesIndices, py::return_value_policy<py::copy_const_reference>())
        .def("allEquilibrium", &Partition::allEquilibrium)
        .def("allKinetic", &Partition::allKinetic)
        .def("allEquilibriumExcept", &Partition::allEquilibriumExcept)
        .def("allKineticExcept", &Partition::allKineticExcept)
        .staticmethod("allEquilibrium")
        .staticmethod("allKinetic")
        .staticmethod("allEquilibriumExcept")
        .staticmethod("allKineticExcept")
		;
}

} // namespace Reaktor
