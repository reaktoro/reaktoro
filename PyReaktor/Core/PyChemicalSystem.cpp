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

#include "PyChemicalSystem.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Interfaces/Gems.hpp>
#include <Reaktor/Interfaces/Phreeqx.hpp>

namespace Reaktor {
namespace {

auto createChemicalSystemGems(const Gems& gems) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(gems);
}

auto createChemicalSystemPhreeqx(const Phreeqx& phreeqx) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(phreeqx);
}

} // namespace

auto export_ChemicalSystem() -> void
{
    py::class_<Connectivity>("Connectivity")
        .def_readwrite("element_to_species", &Connectivity::element_to_species)
        .def_readwrite("species_to_elements", &Connectivity::species_to_elements)
        .def_readwrite("species_to_phase", &Connectivity::species_to_phase)
        .def_readwrite("phase_to_species", &Connectivity::phase_to_species)
        .def_readwrite("element_to_phases", &Connectivity::element_to_phases)
        .def_readwrite("phase_to_elements", &Connectivity::phase_to_elements)
        ;

    py::class_<ChemicalSystem, py::bases<Multiphase>>("ChemicalSystem")
        .def(py::init<>())
        .def(py::init<const Multiphase&>())
        .def("__init__", py::make_constructor(createChemicalSystemGems))
        .def("__init__", py::make_constructor(createChemicalSystemPhreeqx))
        .def("connectivity", &ChemicalSystem::connectivity, py::return_value_policy<py::copy_const_reference>())
        .def(py::self_ns::str(py::self_ns::self));
        ;
}

} // namespace Reaktor
