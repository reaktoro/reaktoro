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

// Reaktoro includes
#include <Reaktoro/Common/TimeUtils.hpp>

namespace Reaktoro {

/// The global auxiliary time variable used by macros tic and toc.
static Time __tic;

/// The global auxiliary time variable used by macros longtic and longtoc.
static Time __longtic;

/// Macro to start timing of a sequence of statements.
#define tic(expr) { __tic = time(); expr; }

/// Macro to start timing of a sequence of statements.
#define longtic(expr) { __longtic = time(); expr; }

/// Macro to get the execution time since last call to tic(expr) macro.
#define toc(res) { res += elapsed(__tic); }

/// Macro to get the execution time since last call to longtic(expr) macro.
#define longtoc(res) { res += elapsed(__longtic); }

/// Macro to measure the elapsed time of an expression execution.
#define timeit(expr, res) { tic(); expr; toc(res); }

/// Macro that disables an expression if profiling is disabled.
#define profiling(expr) expr

} // namespace Reaktoro

#else // NOT REAKTORO_PROFILING

namespace Reaktoro {

/// Macro to start timing of a sequence of statements.
#define tic(expr)

/// Macro to start timing of a sequence of statements.
#define longtic(expr)

/// Macro to get the execution time since last call to tic(expr) macro.
#define toc(res)

/// Macro to get the execution time since last call to longtic(expr) macro.
#define longtoc(res)

/// Macro to measure the elapsed time of an expression execution.
#define timeit(expr, res) expr

/// Macro that disables an expression if profiling is disabled.
#define profiling(expr)

} // namespace Reaktoro

#endif // REAKTORO_PROFILING
