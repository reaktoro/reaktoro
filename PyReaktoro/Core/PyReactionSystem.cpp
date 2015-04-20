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

#include "PyChemicalSystem.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

auto export_ReactionSystem() -> void
{
    py::class_<ReactionSystemModel>("ReactionSystemModel")
        .def_readwrite("lnk", &ReactionSystemModel::lnk)
        .def_readwrite("standard_gibbs_energy", &ReactionSystemModel::standard_gibbs_energy)
        .def_readwrite("standard_helmholtz_energy", &ReactionSystemModel::standard_helmholtz_energy)
        .def_readwrite("standard_internal_energy", &ReactionSystemModel::standard_internal_energy)
        .def_readwrite("standard_enthalpy", &ReactionSystemModel::standard_enthalpy)
        .def_readwrite("standard_entropy", &ReactionSystemModel::standard_entropy)
        .def_readwrite("standard_volume", &ReactionSystemModel::standard_volume)
        .def_readwrite("standard_heat_capacity", &ReactionSystemModel::standard_heat_capacity)
        .def_readwrite("rate", &ReactionSystemModel::rate)
        ;

    using return_const_ref = py::return_value_policy<py::copy_const_reference>;

    auto reaction1 = static_cast<const Reaction&(ReactionSystem::*)(Index) const>(&ReactionSystem::reaction);
    auto reaction2 = static_cast<const Reaction&(ReactionSystem::*)(std::string) const>(&ReactionSystem::reaction);

    py::class_<ReactionSystem>("ReactionSystem")
        .def(py::init<>())
        .def(py::init<const std::vector<Reaction>&>())
        .def(py::init<const std::vector<Reaction>&, const ReactionSystemModel&>())
        .def("numReactions", &ReactionSystem::numReactions)
        .def("indexReaction", &ReactionSystem::indexReaction)
        .def("reactions", &ReactionSystem::reactions, return_const_ref())
        .def("reaction", reaction1, return_const_ref())
        .def("reaction", reaction2, return_const_ref())
        .def("stoichiometricMatrix", &ReactionSystem::stoichiometricMatrix, return_const_ref())
        .def("system", &ReactionSystem::system, return_const_ref())
        .def("lnEquilibriumConstants", &ReactionSystem::lnEquilibriumConstants)
        .def("standardGibbsEnergies", &ReactionSystem::standardGibbsEnergies)
        .def("standardHelmholtzEnergies", &ReactionSystem::standardHelmholtzEnergies)
        .def("standardInternalEnergies", &ReactionSystem::standardInternalEnergies)
        .def("standardEnthalpies", &ReactionSystem::standardEnthalpies)
        .def("standardEntropies", &ReactionSystem::standardEntropies)
        .def("standardVolumes", &ReactionSystem::standardVolumes)
        .def("standardHeatCapacities", &ReactionSystem::standardHeatCapacities)
        .def("rates", &ReactionSystem::rates)
        .def("lnReactionQuotients", &ReactionSystem::lnReactionQuotients)
        ;
}

} // namespace Reaktoro
