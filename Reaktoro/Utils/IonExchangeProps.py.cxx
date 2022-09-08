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
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Utils/IonExchangeProps.hpp>
using namespace Reaktoro;

void exportIonExchangeProps(py::module& m)
{
    py::class_<IonExchangeProps>(m, "IonExchangeProps")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalState&>())
        .def(py::init<const ChemicalProps&>())
        .def("update", py::overload_cast<const ChemicalState&>(&IonExchangeProps::update))
        .def("update", py::overload_cast<const ChemicalProps&>(&IonExchangeProps::update))
        .def("elementAmounts", &IonExchangeProps::elementAmounts)
        .def("elementAmount", &IonExchangeProps::elementAmount)
        .def("speciesAmount", &IonExchangeProps::speciesAmount)
        .def("speciesAmounts", &IonExchangeProps::speciesAmounts)
        .def("speciesEquivalents", &IonExchangeProps::speciesEquivalents)
        .def("speciesEquivalent", &IonExchangeProps::speciesEquivalent)
        .def("speciesEquivalentFractions", &IonExchangeProps::speciesEquivalentFractions)
        .def("speciesEquivalentFraction", &IonExchangeProps::speciesEquivalentFraction)
        .def("speciesActivityCoefficientsLg", &IonExchangeProps::speciesActivityCoefficientsLg)
        .def("speciesActivityCoefficientLg", &IonExchangeProps::speciesActivityCoefficientLg)
        .def("phase", &IonExchangeProps::phase, return_internal_ref)
        .def("output", py::overload_cast<std::ostream&>(&IonExchangeProps::output, py::const_))
        .def("output", py::overload_cast<const String&>(&IonExchangeProps::output, py::const_))
        .def("__repr__", [](const IonExchangeProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
