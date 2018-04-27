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

#include "PyPhreeqcEditor.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Interfaces/Phreeqc.hpp>
#include <Reaktoro/Interfaces/PhreeqcEditor.hpp>

namespace Reaktoro {

auto export_PhreeqcEditor() -> void
{
    py::class_<PhreeqcEditor>("PhreeqcEditor")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("setDatabase", &PhreeqcEditor::setDatabase)
        .def("setAqueousPhase", &PhreeqcEditor::setAqueousPhase)
        .def("setGaseousPhase", &PhreeqcEditor::setGaseousPhase)
        .def("setMineralPhases", &PhreeqcEditor::setMineralPhases)
        ;
}

} // namespace Reaktoro
