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

// PyReaktoro includes
#include <PyReaktoro/Optimization/PyNonlinearOptions.hpp>
#include <PyReaktoro/Optimization/PyOptimumMethod.hpp>
#include <PyReaktoro/Optimization/PyOptimumOptions.hpp>
#include <PyReaktoro/Optimization/PyOptimumResult.hpp>
#include <PyReaktoro/Optimization/PyOptimumState.hpp>

namespace Reaktoro {

inline auto export_Optimization() -> void
{
    export_NonlinearOptions();
    export_OptimumMethod();
    export_OptimumOptions();
    export_OptimumResult();
    export_OptimumState();
}

} // namespace Reaktoro
