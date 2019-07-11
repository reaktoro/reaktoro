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
#include <Reaktoro/Thermodynamics/Species/FluidSpecies.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

namespace Reaktoro {

void exportDatabase(py::module& m)
{
    auto aqueousSpecies1 = static_cast<std::vector<AqueousSpecies>(Database::*)()>(&Database::aqueousSpecies);
    auto aqueousSpecies2 = static_cast<const AqueousSpecies&(Database::*)(std::string) const>(&Database::aqueousSpecies);

    auto fluidSpecies1 = static_cast<std::vector<FluidSpecies>(Database::*)()>(&Database::liquidSpecies);
    auto fluidSpecies2 = static_cast<const FluidSpecies&(Database::*)(std::string) const>(&Database::liquidSpecies);

    auto gaseousSpecies1 = static_cast<std::vector<FluidSpecies>(Database::*)()>(&Database::gaseousSpecies);
    auto gaseousSpecies2 = static_cast<const FluidSpecies&(Database::*)(std::string) const>(&Database::gaseousSpecies);

    auto liquidSpecies1 = static_cast<std::vector<FluidSpecies>(Database::*)()>(&Database::liquidSpecies);
    auto liquidSpecies2 = static_cast<const FluidSpecies&(Database::*)(std::string) const>(&Database::liquidSpecies);

    auto mineralSpecies1 = static_cast<std::vector<MineralSpecies>(Database::*)()>(&Database::mineralSpecies);
    auto mineralSpecies2 = static_cast<const MineralSpecies&(Database::*)(std::string) const>(&Database::mineralSpecies);

    py::class_<Database>(m, "Database")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("elements", &Database::elements)
		.def("addElement", &Database::addElement)
        .def("aqueousSpecies", aqueousSpecies1)
        .def("aqueousSpecies", aqueousSpecies2, py::return_value_policy::reference_internal)
		.def("addAqueousSpecies", &Database::addAqueousSpecies)
        .def("fluidSpecies", fluidSpecies1)
        .def("fluidSpecies", fluidSpecies2)
        .def("addFluidSpecies", &Database::addFluidSpecies)
        .def("gaseousSpecies", gaseousSpecies1)
        .def("gaseousSpecies", gaseousSpecies2, py::return_value_policy::reference_internal)
		.def("addGaseousSpecies", &Database::addGaseousSpecies)
        .def("liquidSpecies", liquidSpecies1)
        .def("liquidSpecies", liquidSpecies2, py::return_value_policy::reference_internal)
        .def("addLiquidSpecies", &Database::addLiquidSpecies)
        .def("mineralSpecies", mineralSpecies1)
        .def("mineralSpecies", mineralSpecies2, py::return_value_policy::reference_internal)
		.def("addMineralSpecies", &Database::addMineralSpecies)
        .def("containsAqueousSpecies", &Database::containsAqueousSpecies)
        .def("containsFluidSpecies", &Database::containsFluidSpecies)
        .def("containsGaseousSpecies", &Database::containsGaseousSpecies)
        .def("containsLiquidSpecies", &Database::containsLiquidSpecies)
        .def("containsMineralSpecies", &Database::containsMineralSpecies)
        .def("aqueousSpeciesWithElements", &Database::aqueousSpeciesWithElements)
        .def("fluidSpeciesWithElements", &Database::fluidSpeciesWithElements)
        .def("gaseousSpeciesWithElements", &Database::gaseousSpeciesWithElements)
        .def("liquidSpeciesWithElements", &Database::liquidSpeciesWithElements)
        .def("mineralSpeciesWithElements", &Database::mineralSpeciesWithElements)
        ;
}

} // namespace Reaktoro
