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

#include "PySpecies.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Species.hpp>

// PyReator includes
#include <PyReaktor/Utils/PyConverters.hpp>

namespace Reaktor {

auto export_Species() -> void
{
    py::class_<Species>("Species")
        .def(py::init<>())
        .def("name", &Species::name, py::return_value_policy<py::copy_const_reference>())
        .def("formula", &Species::formula, py::return_value_policy<py::copy_const_reference>())
        .def("elements", &Species::elements, py::return_value_policy<py::copy_const_reference>())
        .def("elementCoefficient", &Species::elementCoefficient)
        .def("molarMass", &Species::molarMass)
        .def("standardGibbsEnergy", &Species::standardGibbsEnergy)
        .def("standardHelmholtzEnergy", &Species::standardHelmholtzEnergy)
        .def("standardInternalEnergy", &Species::standardInternalEnergy)
        .def("standardEnthalpy", &Species::standardEnthalpy)
        .def("standardEntropy", &Species::standardEntropy)
        .def("standardVolume", &Species::standardVolume)
        .def("standardHeatCapacity", &Species::standardHeatCapacity)
        ;

    export_std_vector<Species>("SpeciesVector");
}

} // namespace Reaktor
