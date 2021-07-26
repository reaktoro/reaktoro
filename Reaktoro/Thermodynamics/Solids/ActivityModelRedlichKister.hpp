// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Return the activity model for a binary solid solution phase based on Redlich-Kister model.
/// The Redlich-Kister model calculates the activity coefficient of the end-members in a
/// solid solution using the equations:
/// @f[\ln\gamma_{1}=x_{2}^{2}[a_{0}+a_{1}(3x_{1}-x_{2})+a_{2}(x_{1}-x_{2})(5x_{1}-x_{2})]@f]
/// and
/// @f[\ln\gamma_{2}=x_{1}^{2}[a_{0}-a_{1}(3x_{2}-x_{1})+a_{2}(x_{2}-x_{1})(5x_{2}-x_{1})]@f].
/// The parameters @f$a_0@f$, @f$a_1@f$, and @f$a_2@f$ must be provided.
/// Set them to zero if not needed.
/// @param mixture The mineral mixture
/// @param a0 The Redlich-Kister parameter a0
/// @param a1 The Redlich-Kister parameter a1
/// @param a2 The Redlich-Kister parameter a2
/// @return The equation of state function for the mineral phase
/// @see ActivityModel
auto ActivityModelRedlichKister(real a0, real a1, real a2) -> ActivityModelGenerator;

} // namespace Reaktoro
