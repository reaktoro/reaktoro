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
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

// PyReator includes
#include <PyReaktoro/Utils/PyConverters.hpp>

namespace Reaktoro {

auto export_ChemicalPlot() -> void
{
    auto y1 = static_cast<void(ChemicalPlot::*)(std::vector<std::string>)>(&ChemicalPlot::y);
    auto y2 = static_cast<void(ChemicalPlot::*)(std::string)>(&ChemicalPlot::y);

    auto legend1 = static_cast<void(ChemicalPlot::*)(std::vector<std::string>)>(&ChemicalPlot::legend);
    auto legend2 = static_cast<void(ChemicalPlot::*)(std::string)>(&ChemicalPlot::legend);

    auto lshift = static_cast<ChemicalPlot&(ChemicalPlot::*)(std::string)>(&ChemicalPlot::operator<<);

    py::class_<ChemicalPlot>("ChemicalPlot")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def("x", &ChemicalPlot::x)
        .def("y", y1)
        .def("y", y2)
        .def("legend", legend1)
        .def("legend", legend2)
        .def("frequency", &ChemicalPlot::frequency)
        .def("__lshift__", lshift, py::return_internal_reference<>())
        .def("open", &ChemicalPlot::open)
        .def("update", &ChemicalPlot::update)
        ;

    export_std_vector<ChemicalPlot>("ChemicalPlotVector");
}

} // namespace Reaktoro
