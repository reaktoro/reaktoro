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

#pragma once

// C++ includes
#include <any>
#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

/// The type used to represent indices and unsigned integers in the library.
using Index = std::size_t;

/// The type that represents a collection of indices.
using Indices = std::vector<Index>;

/// Convenient alias for `std::string`.
using String = std::string;

/// Convenient alias for `std::vector<String>`.
using Strings = std::vector<std::string>;

/// The type used to accept either a name or an index.
using StringOrIndex = std::variant<Index, int, std::string>;

/// Convenient alias for `std::array<T, N>`.
template<typename T, std::size_t N>
using Array = std::array<T, N>;

/// Convenient alias for `std::vector<T>`.
template<typename T>
using Vec = std::vector<T>;

/// Convenient alias for `std::unordered_map<Key, T>`.
template<typename Key, typename T>
using Map = std::unordered_map<Key, T>;

// Conveniet alias to `std::unordered_set<T>`
template<typename T>
using Set = std::unordered_set<T>;

/// Convenient alias for `std::pair<T, U>`.
template<typename T, typename U>
using Pair = std::pair<T, U>;

/// Convenient alias for `std::vector<std::pair<T, U>>`.
template<typename T, typename U>
using Pairs = Vec<Pair<T, U>>;

/// Convenient alias for `std::tuple<Args...>`.
template<typename... Args>
using Tuple = std::tuple<Args...>;

/// Convenient alias for `std::unique_ptr<T>`.
template<typename T>
using Ptr = std::unique_ptr<T>;

/// Convenient alias for `std::shared_ptr<T>`.
template<typename T>
using SharedPtr = std::shared_ptr<T>;

/// Convenient alias for `std::function<R(Args...)>`.
template<typename F>
using Fn = std::function<F>;

/// Convenient alias for `std::optional`.
template<typename T>
using Optional = std::optional<T>;

/// Convenient alias for `std::any`.
using Any = std::any;

} // namespace Reaktoro
