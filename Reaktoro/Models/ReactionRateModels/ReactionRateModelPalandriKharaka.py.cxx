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
        .def_readwrite("formula",  &ReactionRateModelParamsPalandriKharaka::Catalyst::formula , "The chemical formula of the species that participates as a catalyst.")
        .def_readwrite("property", &ReactionRateModelParamsPalandriKharaka::Catalyst::property, "The symbol of the species property that acts as a catalyser. The options are `a` for the activity of the species, and `P` for its partial pressure (in which case the species must be an existing gas in the system). Defaults to 'a'.")
        .def_readwrite("power",    &ReactionRateModelParamsPalandriKharaka::Catalyst::power   , "The power of the property that affects the rate of mineral reaction. Defaults to 0.0.")
        ;

    py::class_<ReactionRateModelParamsPalandriKharaka::Mechanism>(m, "ReactionRateModelParamsPalandriKharakaMechanism")
        .def(py::init<>())
        .def_readwrite("name",      &ReactionRateModelParamsPalandriKharaka::Mechanism::name     , "The classifying name of the mineral reaction mechanism (e.g., `Acid`, `Neutral`, `Base`, `Carbonate`).")
        .def_readwrite("lgk",       &ReactionRateModelParamsPalandriKharaka::Mechanism::lgk      , "The kinetic rate constant of the mineral reaction at 298.15 K (in lg mol/(m2*s)).")
        .def_readwrite("E",         &ReactionRateModelParamsPalandriKharaka::Mechanism::E        , "The Arrhenius activation energy of the mineral reaction (in kJ/mol).")
        .def_readwrite("p",         &ReactionRateModelParamsPalandriKharaka::Mechanism::p        , "The empirical and dimensionless power parameter *p*.")
        .def_readwrite("q",         &ReactionRateModelParamsPalandriKharaka::Mechanism::q        , "The empirical and dimensionless power parameter *q*.")
        .def_readwrite("catalysts", &ReactionRateModelParamsPalandriKharaka::Mechanism::catalysts, "The catalysts of the mineral reaction.")
        ;

    py::class_<ReactionRateModelParamsPalandriKharaka>(m, "ReactionRateModelParamsPalandriKharaka")
        .def(py::init<>())
        .def_readwrite("mineral",    &ReactionRateModelParamsPalandriKharaka::mineral   , "The name of the mineral (e.g., `Dolomite`).")
        .def_readwrite("othernames", &ReactionRateModelParamsPalandriKharaka::othernames, "The alternative names of the mineral (e.g., `Dolomite,ord`, `Dolomite,ordered`).")
        .def_readwrite("mechanisms", &ReactionRateModelParamsPalandriKharaka::mechanisms, "The reaction mechanisms considered in the mineral dissolution/precipitation rate model.")
        ;

    m.def("ReactionRateModelPalandriKharaka", py::overload_cast<>(ReactionRateModelPalandriKharaka));
    m.def("ReactionRateModelPalandriKharaka", py::overload_cast<Params const&>(ReactionRateModelPalandriKharaka));
    m.def("ReactionRateModelPalandriKharaka", py::overload_cast<ReactionRateModelParamsPalandriKharaka const&>(ReactionRateModelPalandriKharaka));
    m.def("ReactionRateModelPalandriKharaka", py::overload_cast<Vec<ReactionRateModelParamsPalandriKharaka> const&>(ReactionRateModelPalandriKharaka));
}
