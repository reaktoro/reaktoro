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
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/ElementalComposition.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Serialization/Core.hpp>
using namespace Reaktoro;

void exportSerializationCore(py::module& m)
{
    py::implicitly_convertible<Data, AggregateState>();
    py::implicitly_convertible<AggregateState, Data>();

    py::implicitly_convertible<Data, ChemicalFormula>();
    py::implicitly_convertible<ChemicalFormula, Data>();

    py::implicitly_convertible<Data, ChemicalSystem>();
    py::implicitly_convertible<ChemicalSystem, Data>();

    py::implicitly_convertible<Data, Database>();
    py::implicitly_convertible<Database, Data>();

    py::implicitly_convertible<Data, Element>();
    py::implicitly_convertible<Element, Data>();

    py::implicitly_convertible<Data, ElementList>();
    py::implicitly_convertible<ElementList, Data>();

    py::implicitly_convertible<Data, ElementalComposition>();
    py::implicitly_convertible<ElementalComposition, Data>();

    py::implicitly_convertible<Data, FormationReaction>();
    py::implicitly_convertible<FormationReaction, Data>();

    py::implicitly_convertible<Data, Param>();
    py::implicitly_convertible<Param, Data>();

    py::implicitly_convertible<Data, Phase>();
    py::implicitly_convertible<Phase, Data>();

    py::implicitly_convertible<Data, ReactionStandardThermoModel>();
    py::implicitly_convertible<ReactionStandardThermoModel, Data>();

    py::implicitly_convertible<Data, Species>();
    py::implicitly_convertible<Species, Data>();

    py::implicitly_convertible<Data, SpeciesList>();
    py::implicitly_convertible<SpeciesList, Data>();

    py::implicitly_convertible<Data, StandardThermoModel>();
    py::implicitly_convertible<StandardThermoModel, Data>();
}
