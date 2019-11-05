// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

void exportChemicalOutput(py::module& m)
{
    auto filename1 = static_cast<void (ChemicalOutput::*)(std::string)>(&ChemicalOutput::filename);
    auto filename2 = static_cast<std::string (ChemicalOutput::*)() const>(&ChemicalOutput::filename);

    auto suffix1 = static_cast<void (ChemicalOutput::*)(std::string)>(&ChemicalOutput::suffix);
    auto suffix2 = static_cast<std::string (ChemicalOutput::*)() const>(&ChemicalOutput::suffix);

    auto add1 = static_cast<void (ChemicalOutput::*)(std::string)>(&ChemicalOutput::add);
    auto add2 = static_cast<void (ChemicalOutput::*)(std::string, std::string)>(&ChemicalOutput::add);

    auto attach1 = static_cast<void (ChemicalOutput::*)(int)>(&ChemicalOutput::attach);
    auto attach2 = static_cast<void (ChemicalOutput::*)(double)>(&ChemicalOutput::attach);
    auto attach3 = static_cast<void (ChemicalOutput::*)(std::string)>(&ChemicalOutput::attach);

    py::class_<ChemicalOutput>(m, "ChemicalOutput")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def("filename", filename1)
        .def("filename", filename2)
        .def("suffix", suffix1)
        .def("suffix", suffix2)
        .def("basename", &ChemicalOutput::basename)
        .def("extension", &ChemicalOutput::extension)
        .def("add", add1)
        .def("add", add2)
        .def("attachments", &ChemicalOutput::attachments)
        .def("attach", attach1)
        .def("attach", attach2)
        .def("attach", attach3)
        .def("precision", &ChemicalOutput::precision)
        .def("scientific", &ChemicalOutput::scientific)
        .def("terminal", &ChemicalOutput::terminal)
        .def("quantities", &ChemicalOutput::quantities)
        .def("headings", &ChemicalOutput::headings)
        .def("open", &ChemicalOutput::open)
        .def("update", &ChemicalOutput::update)
        .def("close", &ChemicalOutput::close);
}

} // namespace Reaktoro
