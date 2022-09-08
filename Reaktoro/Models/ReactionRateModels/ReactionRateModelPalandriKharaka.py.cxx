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
#include <Reaktoro/Models/ReactionRateModels/ReactionRateModelPalandriKharaka.hpp>
using namespace Reaktoro;

void exportReactionRateModelPalandriKharaka(py::module& m)
{
    py::class_<ReactionRateModelParamsPalandriKharaka::Catalyst>(m, "ReactionRateModelParamsPalandriKharakaCatalyst")
        .def(py::init<>())
        .def_readwrite("formula",  &ReactionRateModelParamsPalandriKharaka::Catalyst::formula)
        .def_readwrite("property", &ReactionRateModelParamsPalandriKharaka::Catalyst::property)
        .def_readwrite("power",    &ReactionRateModelParamsPalandriKharaka::Catalyst::power)
        ;

    py::class_<ReactionRateModelParamsPalandriKharaka::Mechanism>(m, "ReactionRateModelParamsPalandriKharakaMechanism")
        .def(py::init<>())
        .def_readwrite("name",      &ReactionRateModelParamsPalandriKharaka::Mechanism::name)
        .def_readwrite("lgk",       &ReactionRateModelParamsPalandriKharaka::Mechanism::lgk)
        .def_readwrite("E",         &ReactionRateModelParamsPalandriKharaka::Mechanism::E)
        .def_readwrite("p",         &ReactionRateModelParamsPalandriKharaka::Mechanism::p)
        .def_readwrite("q",         &ReactionRateModelParamsPalandriKharaka::Mechanism::q)
        .def_readwrite("catalysts", &ReactionRateModelParamsPalandriKharaka::Mechanism::catalysts)
        ;

    py::class_<ReactionRateModelParamsPalandriKharaka>(m, "ReactionRateModelParamsPalandriKharaka")
        .def(py::init<>())
        .def_readwrite("names",      &ReactionRateModelParamsPalandriKharaka::names)
        .def_readwrite("mechanisms", &ReactionRateModelParamsPalandriKharaka::mechanisms)
        ;

    m.def("ReactionRateModelPalandriKharaka", ReactionRateModelPalandriKharaka);
}
