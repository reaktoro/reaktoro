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
#include <Reaktoro/Serialization/Core.YAML.hpp>
using namespace Reaktoro;

void exportSerializationCoreYAML(py::module& m)
{
    py::implicitly_convertible<yaml, AggregateState>();
    py::implicitly_convertible<AggregateState, yaml>();

    py::implicitly_convertible<yaml, ChemicalFormula>();
    py::implicitly_convertible<ChemicalFormula, yaml>();

    py::implicitly_convertible<yaml, ChemicalSystem>();
    py::implicitly_convertible<ChemicalSystem, yaml>();

    py::implicitly_convertible<yaml, Database>();
    py::implicitly_convertible<Database, yaml>();

    py::implicitly_convertible<yaml, Element>();
    py::implicitly_convertible<Element, yaml>();

    py::implicitly_convertible<yaml, ElementList>();
    py::implicitly_convertible<ElementList, yaml>();

    py::implicitly_convertible<yaml, ElementalComposition>();
    py::implicitly_convertible<ElementalComposition, yaml>();

    py::implicitly_convertible<yaml, FormationReaction>();
    py::implicitly_convertible<FormationReaction, yaml>();

    py::implicitly_convertible<yaml, Param>();
    py::implicitly_convertible<Param, yaml>();

    py::implicitly_convertible<yaml, Phase>();
    py::implicitly_convertible<Phase, yaml>();

    py::implicitly_convertible<yaml, ReactionThermoModel>();
    py::implicitly_convertible<ReactionThermoModel, yaml>();

    py::implicitly_convertible<yaml, Species>();
    py::implicitly_convertible<Species, yaml>();

    py::implicitly_convertible<yaml, SpeciesList>();
    py::implicitly_convertible<SpeciesList, yaml>();

    py::implicitly_convertible<yaml, StandardThermoModel>();
    py::implicitly_convertible<StandardThermoModel, yaml>();
}
