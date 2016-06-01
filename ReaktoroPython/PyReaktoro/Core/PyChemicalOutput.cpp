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
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

auto export_ChemicalOutput() -> void
{
    auto addData1 = static_cast<void(ChemicalOutput::*)(std::string)>(&ChemicalOutput::addData);
    auto addData2 = static_cast<void(ChemicalOutput::*)(std::string, std::string)>(&ChemicalOutput::addData);

    py::class_<ChemicalOutput>("ChemicalOutput")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def("setOutputFile", &ChemicalOutput::setOutputFile)
        .def("addData", addData1)
        .def("addData", addData2)
        .def("enableTerminalOutput", &ChemicalOutput::enableTerminalOutput)
        .def("open", &ChemicalOutput::open)
        .def("update", &ChemicalOutput::update)
        .def("open", &ChemicalOutput::close)
        ;
}

} // namespace Reaktoro
