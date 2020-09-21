// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/Database.hpp>
using namespace Reaktoro;

void exportDatabase(py::module& m)
{
    auto addSpecies1 = [](Database& self, const Species& species)
    {
        self.addSpecies(species);
    };

    auto addSpecies2 = [](Database& self, const Vec<Species>& species)
    {
        self.addSpecies(species);
    };

    py::class_<Database>(m, "Database")
        .def(py::init<>())
        .def("clear", &Database::clear)
        .def("addSpecies", addSpecies1)
        .def("addSpecies", addSpecies2)
        .def("attachData", &Database::attachData)
        .def("elements", &Database::elements)
        .def("species", &Database::species)
        .def("speciesWithAggregateState", &Database::speciesWithAggregateState)
        .def("attachedData", &Database::attachedData)
        ;
}
