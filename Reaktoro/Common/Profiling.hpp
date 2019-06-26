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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/TimeUtils.hpp>

#ifdef REAKTORO_PROFILING

#define profiling(expr) expr

#else // NOT REAKTORO_PROFILING

#define profiling(expr) ((void) (0))

#endif // REAKTORO_PROFILING


#ifdef REAKTORO_PROFILING

namespace Reaktoro {

/// A global time variable for profiling codes.
static Time __profiling_timer;

/// Macro to time the execution of an expression.
#define tic(expr) __profiling_timer = time(); expr

/// Get the elapsed time since the last call to timeit(expr) (in units of s)
inline auto toc(double& time) -> void
{
    time = elapsed(__profiling_timer);
}

} // namespace Reaktoro

#else // NOT REAKTORO_PROFILING

namespace Reaktoro {

/// Macro to time the execution of an expression.
#define tic(expr) expr

/// Get the elapsed time since the last call to timeit(expr) (in units of s)
constexpr auto toc(double& time) -> void
{
}

} // namespace Reaktoro

#endif // REAKTORO_PROFILING
