// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <string>
#include <vector>
#include <memory>

namespace Reaktoro {

/// Convenient alias for `std::string`.
using String = std::string;

/// Convenient alias for `std::vector<std::string>`.
using Strings = std::vector<std::string>;

/// Convenient alias for `std::unique_ptr<T>`.
template<typename T>
using Ptr = std::unique_ptr<T>;

/// Convenient alias for `std::shared_ptr<T>`.
template<typename T>
using SharedPtr = std::shared_ptr<T>;

} // namespace Reaktoro
