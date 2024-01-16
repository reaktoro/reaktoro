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
#include <Reaktoro/Common/TimeUtils.hpp>

namespace Reaktoro {

#ifdef REAKTORO_DISABLE_PROFILING

/// Macro to start timing of a sequence of statements.
#define tic(id)

/// Macro to get the execution time since last call to tic(expr) macro.
#define toc(id) 0.0

/// Macro to measure the elapsed time of an expression execution.
#define timeit(expr, res) { expr; res 0.0; }

#else

/// Macro to start timing of a sequence of statements.
#define tic(id) Time __start_time ## id = time();

/// Macro to get the execution time since last call to tic(expr) macro.
#define toc(id) elapsed(__start_time ## id)

/// Macro to measure the elapsed time of an expression execution.
#define timeit(expr, res) { tic(__##__LINE__); expr; res elapsed(__start_time##__##__LINE__); }

#endif // REAKTORO_DISABLE_PROFILING

} // namespace Reaktoro
