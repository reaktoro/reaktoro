// Reaktoro is a unified framework for modeling chemically reactive systems.
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
#include <Reaktoro/Common/Json.hpp>

namespace Reaktoro {

// Forward declarations (struct)
struct KineticResult;
struct KineticTiming;
struct SmartKineticResult;
struct SmartKineticTiming;

// Json converters for KineticResult
void to_json(json& j, const KineticResult& obj);
void from_json(const json& j, KineticResult& obj);

// Json converters for KineticTiming
void to_json(json& j, const KineticTiming& obj);
void from_json(const json& j, KineticTiming& obj);

// Json converters for SmartKineticResult
void to_json(json& j, const SmartKineticResult& obj);
void from_json(const json& j, SmartKineticResult& obj);

// Json converters for SmartKineticTiming
void to_json(json& j, const SmartKineticTiming& obj);
void from_json(const json& j, SmartKineticTiming& obj);

} // namespace Reaktoro
