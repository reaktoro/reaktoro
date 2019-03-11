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
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

namespace Reaktoro {

void exportDatabase(py::module& m)
{
    auto aqueousSpecies1 = static_cast<std::vector<AqueousSpecies>(Database::*)()>(&Database::aqueousSpecies);
    auto aqueousSpecies2 = static_cast<const AqueousSpecies&(Database::*)(std::string) const>(&Database::aqueousSpecies);

    auto gaseousSpecies1 = static_cast<std::vector<GaseousSpecies>(Database::*)()>(&Database::gaseousSpecies);
    auto gaseousSpecies2 = static_cast<const GaseousSpecies&(Database::*)(std::string) const>(&Database::gaseousSpecies);

    auto mineralSpecies1 = static_cast<std::vector<MineralSpecies>(Database::*)()>(&Database::mineralSpecies);
    auto mineralSpecies2 = static_cast<const MineralSpecies&(Database::*)(std::string) const>(&Database::mineralSpecies);

    py::class_<Database>(m, "Database")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("elements", &Database::elements)
        .def("aqueousSpecies", aqueousSpecies1)
        .def("aqueousSpecies", aqueousSpecies2, py::return_value_policy::reference_internal)
        .def("gaseousSpecies", gaseousSpecies1)
        .def("gaseousSpecies", gaseousSpecies2, py::return_value_policy::reference_internal)
        .def("mineralSpecies", mineralSpecies1)
        .def("mineralSpecies", mineralSpecies2, py::return_value_policy::reference_internal)
        .def("containsAqueousSpecies", &Database::containsAqueousSpecies)
        .def("containsGaseousSpecies", &Database::containsGaseousSpecies)
        .def("containsMineralSpecies", &Database::containsMineralSpecies)
        .def("aqueousSpeciesWithElements", &Database::aqueousSpeciesWithElements)
        .def("gaseousSpeciesWithElements", &Database::gaseousSpeciesWithElements)
        .def("mineralSpeciesWithElements", &Database::mineralSpeciesWithElements)
        ;
}

} // namespace Reaktoro
