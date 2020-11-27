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
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/SmartKineticOptions.hpp>

namespace Reaktoro {

/// The structure with options for the kinetic path calculations.
struct KineticPathOptions
{
    /// The options for the equilibrium calculations.
    EquilibriumOptions equilibrium;

    /// The options for the smart equilibrium calculations.
    SmartEquilibriumOptions smart_equilibrium;

    /// The options for the kinetic calculations.
    KineticOptions kinetics;

    /// The options for the smart kinetic calculations.
    SmartKineticOptions smart_kinetics;

    /// The boolean flag that indicates whether smart equilibrium solver should be used.
    bool use_smart_equilibrium_solver = false;

    /// The boolean flag that indicates whether smart kinetic solver should be used.
    bool use_smart_kinetic_solver = false;
};

} // namespace Reaktoro
