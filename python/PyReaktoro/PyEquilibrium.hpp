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
#include <PyReaktoro/Equilibrium/PyEquilibriumCompositionProblem.hpp>
#include <PyReaktoro/Equilibrium/PyEquilibriumInverseProblem.hpp>
#include <PyReaktoro/Equilibrium/PyEquilibriumOptions.hpp>
#include <PyReaktoro/Equilibrium/PyEquilibriumPath.hpp>
#include <PyReaktoro/Equilibrium/PyEquilibriumProblem.hpp>
#include <PyReaktoro/Equilibrium/PyEquilibriumResult.hpp>
#include <PyReaktoro/Equilibrium/PyEquilibriumSensitivity.hpp>
#include <PyReaktoro/Equilibrium/PyEquilibriumSolver.hpp>
#include <PyReaktoro/Equilibrium/PyEquilibriumUtils.hpp>

namespace Reaktoro {

inline auto export_Equilibrium() -> void
{
    export_EquilibriumOptions();
    export_EquilibriumPath();
    export_EquilibriumCompositionProblem();
    export_EquilibriumInverseProblem();
    export_EquilibriumProblem();
    export_EquilibriumResult();
    export_EquilibriumSensitivity();
    export_EquilibriumSolver();
    export_EquilibriumUtils();
}

} // namespace Reaktoro
