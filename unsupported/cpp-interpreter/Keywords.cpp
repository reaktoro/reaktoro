// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <unsupported/cpp-interpreter/Keywords.hpp>

namespace Reaktoro {
namespace kwd {

ValueUnits::ValueUnits()
{}

ValueUnits::ValueUnits(double value, std::string units)
: value(value), units(units)
{}

ValueUnits::ValueUnits(std::string str)
{
    auto words = split(str);
    Assert(words.size() == 2, "Could not create a ValueUnits instance from `" + str + "`.",
        "Expecting two words in the format `value units`, e.g., `300 kelvin`, `50 moles`");
    value = tofloat(words[0]);
    units = words[1];
}

EntityValue::EntityValue()
{}

EntityValue::EntityValue(std::string entity, double value, std::string units)
: entity(entity), value(value)
{}

EntityValue::EntityValue(std::string str)
{
    auto words = split(str);
    Assert(words.size() == 2, "Could not create an EntityValue instance from `" + str + "`.",
        "Expecting two words in the format `entity value`, e.g., `Na 20`");
    entity = words[0];
    value = tofloat(words[1]);
}

EntityValueUnits::EntityValueUnits()
{}

EntityValueUnits::EntityValueUnits(std::string entity, double value, std::string units)
: ValueUnits(value, units), entity(entity)
{}

EntityValueUnits::EntityValueUnits(std::string str)
{
    auto words = split(str);
    Assert(words.size() == 3, "Could not create an EntityValueUnits instance from `" + str + "`.",
        "Expecting three words in the format `entity value units`, e.g., `Calcite 100 g`");
    entity = words[0];
    value = tofloat(words[1]);
    units = words[2];
}

ValueUnitsEntity::ValueUnitsEntity()
{}

ValueUnitsEntity::ValueUnitsEntity(double value, std::string units, std::string entity)
: EntityValueUnits(entity, value, units)
{}

ValueUnitsEntity::ValueUnitsEntity(std::string str)
{
    auto words = split(str);
    Assert(words.size() == 3, "Could not parse `" + str + "`.",
        "Expecting three words in the format `value units entity`, e.g., `1 kg H2O`");
    value = tofloat(words[0]);
    units = words[1];
    entity = words[2];
}

} // namespace kwd
} // namespace Reaktoro
