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

// autodiff includes
#include <autodiff/forward/real.hpp>

namespace Reaktoro {

/// The number type used throughout the library.
using real = autodiff::real;

/// The method that extracts derivatives from real numbers and vectors/arrays of such numbers.
using autodiff::grad;

/// Convenient method to specify which variable autodiff is to be carried out with respect to.
template<typename T>
auto wrt(T& x) -> T& { return x; }

} // namespace Reaktoro
