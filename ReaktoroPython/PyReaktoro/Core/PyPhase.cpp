// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "PyPhase.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>

// PyReator includes
#include <PyReaktoro/Utils/PyConverters.hpp>

namespace Reaktoro {

auto export_Phase() -> void
{
    auto elements1 = static_cast<const std::vector<Element>&(Phase::*)() const>(&Phase::elements);
    auto elements2 = static_cast<std::vector<Element>&(Phase::*)()>(&Phase::elements);

    auto species1 = static_cast<const std::vector<Species>&(Phase::*)() const>(&Phase::species);
    auto species2 = static_cast<std::vector<Species>&(Phase::*)()>(&Phase::species);
    auto species3 = static_cast<const Species&(Phase::*)(Index) const>(&Phase::species);

    auto properties1 = static_cast<ThermoProperties(Phase::*)(double,double) const>(&Phase::properties);
    auto properties2 = static_cast<PhaseChemicalProperties(Phase::*)(double,double,const Vector&) const>(&Phase::properties);

    py::class_<Phase>("Phase")
        .def(py::init<>())
        .def("setName", &Phase::setName)
        .def("setSpecies", &Phase::setSpecies)
        .def("setThermoModel", &Phase::setThermoModel)
        .def("setChemicalModel", &Phase::setChemicalModel)
        .def("setFluid", &Phase::setFluid)
        .def("setSolid", &Phase::setSolid)
        .def("numElements", &Phase::numElements)
        .def("numSpecies", &Phase::numSpecies)
        .def("name", &Phase::name)
        .def("elements", elements1, py::return_internal_reference<>())
        .def("elements", elements2, py::return_internal_reference<>())
        .def("species", species1, py::return_internal_reference<>())
        .def("species", species2, py::return_internal_reference<>())
        .def("species", species3, py::return_internal_reference<>())
        .def("fluid", &Phase::fluid)
        .def("solid", &Phase::solid)
        .def("indexSpecies", &Phase::indexSpecies)
        .def("properties", properties1)
        .def("properties", properties2)
        ;

    export_std_vector<Phase>("PhaseVector");
}

} // namespace Reaktoro
