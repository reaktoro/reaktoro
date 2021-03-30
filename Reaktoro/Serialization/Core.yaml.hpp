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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/AggregateState.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalFormula;
class ChemicalSystem;
class Element;
class ElementalComposition;
class FormationReaction;
class Param;
class Params;
class Phase;
class Species;
template<typename Signature> class Model;

REAKTORO_YAML_ENCODE_DECLARE(AggregateState);
REAKTORO_YAML_DECODE_DECLARE(AggregateState);

REAKTORO_YAML_ENCODE_DECLARE(ChemicalFormula);
REAKTORO_YAML_DECODE_DECLARE(ChemicalFormula);

REAKTORO_YAML_ENCODE_DECLARE(ChemicalSystem);
REAKTORO_YAML_DECODE_DECLARE(ChemicalSystem);

REAKTORO_YAML_ENCODE_DECLARE(Element);
REAKTORO_YAML_DECODE_DECLARE(Element);

REAKTORO_YAML_ENCODE_DECLARE(ElementalComposition);
REAKTORO_YAML_DECODE_DECLARE(ElementalComposition);

REAKTORO_YAML_ENCODE_DECLARE(FormationReaction);
REAKTORO_YAML_DECODE_DECLARE(FormationReaction);

REAKTORO_YAML_ENCODE_DECLARE(Param);
REAKTORO_YAML_DECODE_DECLARE(Param);

REAKTORO_YAML_ENCODE_DECLARE(Params);
REAKTORO_YAML_DECODE_DECLARE(Params);

REAKTORO_YAML_ENCODE_DECLARE(Phase);
REAKTORO_YAML_DECODE_DECLARE(Phase);

REAKTORO_YAML_ENCODE_DECLARE(Species);
REAKTORO_YAML_DECODE_DECLARE(Species);

} // namespace Reaktoro

// REAKTORO_DECLARE_CONVERT_YAML(ChemicalFormula);

// namespace YAML {

// template<>
// struct convert<ChemicalFormula>
// {
//     static auto encode(const ChemicalFormula& obj);
//     static auto decode(const Node& node, ChemicalFormula& obj);
// };

// } // namespace YAML
