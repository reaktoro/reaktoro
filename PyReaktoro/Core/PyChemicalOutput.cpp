// Reaktoro is a C++ library for computational reaction modelling.
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

#include "PyChemicalOutput.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

auto export_ChemicalOutput() -> void
{
    auto data1 = static_cast<void(ChemicalOutput::*)(std::vector<std::string>)>(&ChemicalOutput::data);
    auto data2 = static_cast<void(ChemicalOutput::*)(std::string)>(&ChemicalOutput::data);

    auto header1 = static_cast<void(ChemicalOutput::*)(std::vector<std::string>)>(&ChemicalOutput::header);
    auto header2 = static_cast<void(ChemicalOutput::*)(std::string)>(&ChemicalOutput::header);

    py::class_<ChemicalOutput>("ChemicalOutput")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def("file", &ChemicalOutput::file)
        .def("terminal", &ChemicalOutput::terminal)
        .def("data", data1)
        .def("data", data2)
        .def("header", header1)
        .def("header", header2)
        .def("open", &ChemicalOutput::open)
        .def("update", &ChemicalOutput::update)
        .def("open", &ChemicalOutput::close)
        ;
}

} // namespace Reaktoro
