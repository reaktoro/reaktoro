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

#include "PyChemicalOutput.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

auto export_ChemicalOutput() -> void
{
    py::class_<ChemicalOutput>("ChemicalOutput")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def("file", &ChemicalOutput::file)
        .def("data", &ChemicalOutput::data)
        .def("headings", &ChemicalOutput::headings)
        .def("terminal", &ChemicalOutput::terminal)
        .def("open", &ChemicalOutput::open)
        .def("update", &ChemicalOutput::update)
        .def("open", &ChemicalOutput::close)
        ;

    py::implicitly_convertible<ChemicalOutput, bool>();
}

} // namespace Reaktoro
