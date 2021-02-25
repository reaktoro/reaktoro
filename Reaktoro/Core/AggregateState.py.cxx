// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/AggregateState.hpp>
using namespace Reaktoro;

void exportAggregateState(py::module& m)
{
    py::enum_<AggregateState>(m, "AggregateState")
        .value("Gas", AggregateState::Gas)
        .value("Liquid", AggregateState::Liquid)
        .value("Solid", AggregateState::Solid)
        .value("Plasma", AggregateState::Plasma)
        .value("CondensedPhase", AggregateState::CondensedPhase)
        .value("Fluid", AggregateState::Fluid)
        .value("LiquidCrystal", AggregateState::LiquidCrystal)
        .value("CrystallineSolid", AggregateState::CrystallineSolid)
        .value("AmorphousSolid", AggregateState::AmorphousSolid)
        .value("Vitreous", AggregateState::Vitreous)
        .value("Adsorbed", AggregateState::Adsorbed)
        .value("Monomeric", AggregateState::Monomeric)
        .value("Polymeric", AggregateState::Polymeric)
        .value("SolidSolution", AggregateState::SolidSolution)
        .value("IonExchange", AggregateState::IonExchange)
        .value("Aqueous", AggregateState::Aqueous)
        .value("Undefined", AggregateState::Undefined)
        ;

    m.def("parseAggregateState", parseAggregateState);
    m.def("identifyAggregateState", identifyAggregateState);
}
