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
#include <Reaktoro/Singletons/CriticalProps.hpp>
using namespace Reaktoro;

void exportCriticalProps(py::module& m)
{
    py::class_<SubstanceCriticalPropsData>(m, "SubstanceCriticalPropsData")
        .def(py::init<>())
        .def_readwrite("Tcr", &SubstanceCriticalPropsData::Tcr)
        .def_readwrite("Pcr", &SubstanceCriticalPropsData::Pcr)
        .def_readwrite("omega", &SubstanceCriticalPropsData::omega)
        ;

    py::class_<SubstanceCriticalProps>(m, "SubstanceCriticalProps")
        .def(py::init<const StringList&>())
        .def(py::init<const SubstanceCriticalPropsData&, const StringList&>())
        .def("setTemperature", &SubstanceCriticalProps::setTemperature)
        .def("setPressure", &SubstanceCriticalProps::setPressure)
        .def("setAcentricFactor", &SubstanceCriticalProps::setAcentricFactor)
        .def("names", &SubstanceCriticalProps::names, return_internal_ref)
        .def("temperature", &SubstanceCriticalProps::temperature)
        .def("pressure", &SubstanceCriticalProps::pressure)
        .def("acentricFactor", &SubstanceCriticalProps::acentricFactor)
        .def("data", &SubstanceCriticalProps::data, return_internal_ref)
        ;

    py::class_<CriticalProps, std::unique_ptr<CriticalProps, py::nodelete>>(m, "CriticalProps")
        .def(py::init([]() { return std::unique_ptr<CriticalProps, py::nodelete>(&CriticalProps::instance()); }))
        .def_static("instance", &CriticalProps::instance, py::return_value_policy::reference)
        .def_static("data", &CriticalProps::data, py::return_value_policy::reference)
        .def_static("defaultCriticalProps", &CriticalProps::defaultCriticalProps, py::return_value_policy::reference)
        .def_static("append", &CriticalProps::append)
        .def_static("overwrite", &CriticalProps::overwrite)
        .def_static("setMissingAs", &CriticalProps::setMissingAs)
        .def_static("size", &CriticalProps::size)
        .def_static("find", &CriticalProps::find)
        .def_static("get", py::overload_cast<const String&>(&CriticalProps::get))
        .def_static("get", py::overload_cast<const StringList&>(&CriticalProps::get))
        .def("__getitem__", [](const CriticalProps& self, Index i) { return self.data()[i]; }, py::return_value_policy::reference)
        .def("__iter__", [](const CriticalProps& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()) // keep object alive while iterator exists;
        ;
}
