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

#ifdef REAKTORO_PROFILING

// C++ includes
#include <stack>

// Reaktoro includes
#include <Reaktoro/Common/TimeUtils.hpp>

namespace Reaktoro {

/// The global stack of collected stated time using macro tic(expr).
static std::stack<Time> __tics;

/// The global auxiliary time variable used in macro toc().
static Time __tic;

/// Provide operator>> to transfer toc() result to a double time variable.
struct __ElapsedTime
{
	double secs;
	inline auto operator>>(double& out) const { out = secs; }
    inline operator double() const { return secs; }
	inline auto accumulate(double& out) const { out += secs; }
};

/// Macro to start the timing of an expression execution.
#define tic(expr) __tics.push( time() ); expr;

/// Macro to get the execution time since last call to tic(expr) macro.
#define toc() __tic = __tics.top(); __tics.pop(); __ElapsedTime{elapsed(__tic)}

/// Macro to measure the elapsed time of an expression execution.
#define timeit(expr) tic(); expr; toc()

/// Macro that disables an expression if profiling is disabled.
#define ifprofiling(expr) expr

} // namespace Reaktoro

#else // NOT REAKTORO_PROFILING

namespace Reaktoro {

/// Provide operator>> to transfer toc() result to a double time variable.
struct __ElapsedTime
{
	inline auto operator>>(double& out) const { out = 0.0; }
    inline operator double() const { return 0.0; }
    inline auto accumulate(double& out) const { }
};

/// Macro to start the timing of an expression execution.
#define tic(expr) ((void) (0))

/// Macro to get the execution time since last call to tic(expr) macro.
#define toc() __ElapsedTime{}

/// Macro to measure the elapsed time of an expression execution.
#define timeit(expr) expr; toc()

/// Macro that disables an expression if profiling is disabled.
#define ifprofiling(expr) ((void) (0))

} // namespace Reaktoro

#endif // REAKTORO_PROFILING
