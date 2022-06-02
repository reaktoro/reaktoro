// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Thermodynamics/Surface/DoubleLayerProps.hpp>

using namespace Reaktoro;

void exportDoubleLayerProps(py::module& m)
{
    py::class_<DoubleLayerProps>(m, "DoubleLayerProps")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalState&>())
        .def(py::init<const ChemicalProps&>())
        .def("update", py::overload_cast<const ChemicalState&>(&DoubleLayerProps::update))
        .def("update", py::overload_cast<const ChemicalProps&>(&DoubleLayerProps::update))
        .def("speciesAmount", &DoubleLayerProps::speciesAmount)
        .def("speciesAmounts", &DoubleLayerProps::speciesAmounts)
        .def("speciesMoleFractions", &DoubleLayerProps::speciesMoleFractions)
        .def("speciesMoleFraction", &DoubleLayerProps::speciesMoleFraction)
        .def("output", py::overload_cast<std::ostream&>(&DoubleLayerProps::output, py::const_))
        .def("output", py::overload_cast<const String&>(&DoubleLayerProps::output, py::const_))
        .def("__repr__", [](const DoubleLayerProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
