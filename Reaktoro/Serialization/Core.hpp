// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/Data.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalFormula;
class ChemicalSystem;
class Database;
class Element;
class ElementalComposition;
class FormationReaction;
class Phase;
class Species;
template<typename Signature> class Model;

REAKTORO_DATA_ENCODE_DECLARE(AggregateState);
REAKTORO_DATA_DECODE_DECLARE(AggregateState);

REAKTORO_DATA_ENCODE_DECLARE(ChemicalFormula);
REAKTORO_DATA_DECODE_DECLARE(ChemicalFormula);

REAKTORO_DATA_ENCODE_DECLARE(ChemicalSystem);
REAKTORO_DATA_DECODE_DECLARE(ChemicalSystem);

REAKTORO_DATA_ENCODE_DECLARE(Database);
REAKTORO_DATA_DECODE_DECLARE(Database);

REAKTORO_DATA_ENCODE_DECLARE(Element);
REAKTORO_DATA_DECODE_DECLARE(Element);

REAKTORO_DATA_ENCODE_DECLARE(ElementList);
REAKTORO_DATA_DECODE_DECLARE(ElementList);

REAKTORO_DATA_ENCODE_DECLARE(ElementalComposition);
REAKTORO_DATA_DECODE_DECLARE(ElementalComposition);

REAKTORO_DATA_ENCODE_DECLARE(FormationReaction);
REAKTORO_DATA_DECODE_DECLARE(FormationReaction);

REAKTORO_DATA_ENCODE_DECLARE(Phase);
REAKTORO_DATA_DECODE_DECLARE(Phase);

REAKTORO_DATA_ENCODE_DECLARE(ReactionStandardThermoModel);
REAKTORO_DATA_DECODE_DECLARE(ReactionStandardThermoModel);

REAKTORO_DATA_ENCODE_DECLARE(Species);
REAKTORO_DATA_DECODE_DECLARE(Species);

REAKTORO_DATA_ENCODE_DECLARE(SpeciesList);
REAKTORO_DATA_DECODE_DECLARE(SpeciesList);

REAKTORO_DATA_ENCODE_DECLARE(StandardThermoModel);
REAKTORO_DATA_DECODE_DECLARE(StandardThermoModel);

//--------------------------------------------------------------------------------
// ATTENTION!
//--------------------------------------------------------------------------------
// Linking error if using g++ 12.1 if these encode/decode definitions for
// AggregateState are implemented in cpp file and not here in this header file!
//--------------------------------------------------------------------------------

inline REAKTORO_DATA_ENCODE_DEFINE(AggregateState)
{
    std::stringstream ss;
    ss << obj;
    data = ss.str();
}

inline REAKTORO_DATA_DECODE_DEFINE(AggregateState)
{
    obj = parseAggregateState(data.asString());
}

} // namespace Reaktoro
