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

/// Return the activity model for the ion exchange.
/// @ingroup ActivityModels
auto ActivityModelIonExchange() -> ActivityModelGenerator;

/// Return the Gaines-Thomas activity model for ion exchange.
/// @ingroup ActivityModels
auto ActivityModelIonExchangeGainesThomas() -> ActivityModelGenerator;

//=====================================================================================================================
/// @page ActivityModelIonExchangeGainesThomas Gaines--Thomas activity model
/// The Gaines--Thomas activity model for ion exchange compositions.
/// An instance of this class can be used to control how activities of the
/// ion exchange species are calculated using the Gaines--Thomas activity model.
///
/// The activity of \bold{ion exchange species} are calculated using the
/// equivalence fractions:
///
/// \eqc{\beta_I=-\dfrac{x_{I}Ze_{I}}{\sum_{j} x_{j}Ze_{j},}
///
/// where @eq{x_{I}} and @eq{x_{i}} are species fractions and
/// @eq{Ze_{I}} and @eq{Ze_{i}} are exchanger equivalences (or cation charges)
/// in ion exchange species.

} // namespace Reaktoro
