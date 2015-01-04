// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktor includes
#include <Reaktor/Optimization/OptimumOptions.hpp>

namespace Reaktor {

/// The options for the equilibrium calculations
struct EquilibriumOptions
{
    /// A type to describe the computing options at the end of the equilibrium calculation
    struct Compute
    {
        /// Indicate if the partial derivatives of the equilibrium molar
        /// abundance of the equilibrium species @f$ n @f$ w.r.t.
        /// temperature @f$ T @f$ is to be computed.
        /// By setting this flag to `true`, the partial derivatives
        /// @f$\left.\frac{\partial n}{\partial T}\right|_{P,b}@f$ will be
        /// computed at the end of the equilibrium calculation.
        bool dndt = false;

        /// Indicate if the partial derivatives of the equilibrium molar
        /// abundance of the equilibrium species @f$ n @f$ w.r.t.
        /// pressure @f$ P @f$ is to be computed.
        /// By setting this flag to `true`, the partial derivatives
        /// @f$\left.\frac{\partial n}{\partial P}\right|_{T,b}@f$ will be
        /// computed at the end of the equilibrium calculation.
        bool dndp = false;

        /// Indicate if the partial derivatives of the equilibrium molar
        /// abundance of the equilibrium species @f$ n @f$ w.r.t.
        /// the molar abundance of the elements @f$ b @f$ is to be computed.
        /// By setting this flag to `true`, the partial derivatives
        /// @f$\left.\frac{\partial n}{\partial b}\right|_{T,P}@f$ will be
        /// computed at the end of the equilibrium calculation.
        bool dndb = false;
    };

    /// Construct a default EquilibriumOptions instance
    EquilibriumOptions();

    /// The computing options at the end of the equilibrium calculation
    Compute compute;

    /// The options for the optimisation calculation.
    OptimumOptions optimum;
};

} // namespace Reaktor
