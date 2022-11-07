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
#include <Reaktoro/Thermodynamics/Surface/SurfaceSite.hpp>
#include <Reaktoro/Thermodynamics/Surface/SurfaceSiteProps.hpp>

using namespace Reaktoro;

void exportSurfaceSiteProps(py::module& m)
{
    py::class_<SurfaceSiteProps>(m, "SurfaceSiteProps")
        .def(py::init<const SurfaceSite&, const ChemicalSystem&>())
        .def(py::init<const SurfaceSite&, const ChemicalState&>())
        .def(py::init<const SurfaceSite&, const ChemicalProps&>())
        .def("update", py::overload_cast<const ChemicalState&>(&SurfaceSiteProps::update))
        .def("update", py::overload_cast<const ChemicalProps&>(&SurfaceSiteProps::update))
        .def("elementAmounts", &SurfaceSiteProps::elementAmounts)
        .def("elementAmount", &SurfaceSiteProps::elementAmount)
        .def("speciesAmount", &SurfaceSiteProps::speciesAmount)
        .def("speciesAmounts", &SurfaceSiteProps::speciesAmounts)
        .def("speciesFractions", &SurfaceSiteProps::speciesFractions)
        .def("speciesFraction", &SurfaceSiteProps::speciesFraction)
        .def("charge", &SurfaceSiteProps::charge)
        .def("sigma", &SurfaceSiteProps::sigma)
        .def("output", py::overload_cast<std::ostream&>(&SurfaceSiteProps::output, py::const_))
        .def("output", py::overload_cast<const String&>(&SurfaceSiteProps::output, py::const_))
        .def("__repr__", [](const SurfaceSiteProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
