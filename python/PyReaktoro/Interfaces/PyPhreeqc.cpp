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
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Interfaces/Interface.hpp>
#include <Reaktoro/Interfaces/Phreeqc.hpp>

namespace Reaktoro {

void exportPhreeqc(py::module& m)
{
	auto execute1 = static_cast<void(Phreeqc::*)(std::string,std::string)>(&Phreeqc::execute);
	auto execute2 = static_cast<void(Phreeqc::*)(std::string)>(&Phreeqc::execute);

    py::class_<Phreeqc, Interface>(m, "Phreeqc")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("load", &Phreeqc::load)
        .def("execute", execute1)
        .def("execute", execute2)
        .def("reset", &Phreeqc::reset)
        ;
}

} // namespace Reaktoro
