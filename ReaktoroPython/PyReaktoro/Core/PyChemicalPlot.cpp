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
#include <PyReaktoro/Common/PyConverters.hpp>

namespace Reaktoro {

auto export_ChemicalPlot() -> void
{
    auto addYData1 = static_cast<void(ChemicalPlot::*)(std::string)>(&ChemicalPlot::addYData);
    auto addYData2 = static_cast<void(ChemicalPlot::*)(std::string, std::string)>(&ChemicalPlot::addYData);

    auto lshift = static_cast<ChemicalPlot&(ChemicalPlot::*)(std::string)>(&ChemicalPlot::operator<<);

    py::class_<ChemicalPlot>("ChemicalPlot")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def("setName", &ChemicalPlot::setName)
        .def("setXData", &ChemicalPlot::setXData)
        .def("addYData", addYData1)
        .def("addYData", addYData2)
        .def("setXLabel", &ChemicalPlot::setXLabel)
        .def("setYLabel", &ChemicalPlot::setYLabel)
        .def("setXTics", &ChemicalPlot::setXTics)
        .def("setYTics", &ChemicalPlot::setYTics)
        .def("setXFormat", &ChemicalPlot::setXFormat)
        .def("setYFormat", &ChemicalPlot::setYFormat)
        .def("setXLogscale", &ChemicalPlot::setXLogscale)
        .def("setYLogscale", &ChemicalPlot::setYLogscale)
        .def("setKey", &ChemicalPlot::setKey)
        .def("setRefreshRate", &ChemicalPlot::setRefreshRate)
        .def("enableLegend", &ChemicalPlot::enableLegend)
        .def("__lshift__", lshift, py::return_internal_reference<>())
        .def("open", &ChemicalPlot::open)
        .def("update", &ChemicalPlot::update)
        ;

    export_std_vector<ChemicalPlot>("ChemicalPlotVector");
}

} // namespace Reaktoro
