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
struct ReactiveTransportAnalysis;
struct TransportTiming;
struct TransportResult;

// Json converters for TransportTiming
void to_json(json& j, const TransportTiming& obj);
void from_json(const json& j, TransportTiming& obj);

// Json converters for TransportResult
void to_json(json& j, const TransportResult& obj);
void from_json(const json& j, TransportResult& obj);

// Json converters for ReactiveTransportAnalysis
void to_json(json& j, const ReactiveTransportAnalysis& obj);
void from_json(const json& j, ReactiveTransportAnalysis& obj);

} // namespace Reaktoro
