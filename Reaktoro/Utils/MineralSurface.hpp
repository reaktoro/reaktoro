// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Surfaces.hpp>

namespace Reaktoro {

/// Used to define mineral surfaces that react with aqueous phases.
class MineralSurface : public GeneralSurface
{
public:
    /// Construct a MineralSurface object.
    /// @param mineral The name of the mineral phase in the chemical system.
    explicit MineralSurface(String const& mineral);

    /// Construct a MineralSurface object with a constant or linear surface area model.
    /// If `unitA` is convertible to `m2`, then a constant surface area model is
    /// assigned to the mineral surface (@see SurfaceAreaModelConstant). If
    /// `unitA` is convertible to `m2/mol`, `m2/kg`, or `m2/m3`, then a linear
    /// surface area model is assigned instead (@see SurfaceAreaModelLinear).
    /// Otherwise, an error occurs due to non-conforming unit value for `unitA`.
    /// @param mineral The name of the mineral phase in the chemical system.
    /// @param A The actual or normalized surface area value of the mineral phase.
    /// @param unitA The unit of `A` (must be convertible to `m2`, `m2/mol`, `m2/kg`, or `m2/m3`).
    MineralSurface(String const& mineral, real A, Chars unitA);

    /// Construct a MineralSurface object with a power surface area model.
    /// For more details about the power law used for the surface area model,
    /// see @ref SurfaceAreaModelPower.
    /// @param mineral The name of the mineral phase in the chemical system.
    /// @param A0 The initial surface area of the mineral phase.
    /// @param unitA0 The unit of `A0` (must be convertible to `m2`).
    /// @param q0 The initial quantity of the mineral phase.
    /// @param unitq0 The unit of `q0` (must be convertible to `mol`, `kg`, or `m3`).
    /// @param p The power parameter in the model.
    MineralSurface(String const& mineral, real A0, Chars unitA0, real q0, Chars unitq0, real p);
};

} // namespace Reaktoro
