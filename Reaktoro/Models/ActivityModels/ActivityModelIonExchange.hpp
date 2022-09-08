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
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

/// Return the activity model for the ion exchange.
/// @see @ref PageActivityModelIonExchangeGainesThomas
/// @ingroup ActivityModels
auto ActivityModelIonExchange() -> ActivityModelGenerator;

/// Return the Gaines-Thomas activity model for ion exchange.
/// @see @ref PageActivityModelIonExchangeGainesThomas
/// @ingroup ActivityModels
auto ActivityModelIonExchangeGainesThomas() -> ActivityModelGenerator;

/// Return the Vanselov activity model for ion exchange.
/// @see @ref PageActivityModelIonExchangeVanselov
/// @ingroup ActivityModels
auto ActivityModelIonExchangeVanselow() -> ActivityModelGenerator;

//=====================================================================================================================
/// @page PageActivityModelIonExchangeGainesThomas Gaines--Thomas ion exchange activity model
///
/// The Gaines--Thomas activity model for ion exchange compositions. An
/// instance of this class can be used to control how activities of the ion
/// exchange species are calculated using the Gaines--Thomas activity model.
///
/// The activity of **ion exchange species** are calculated using the
/// either equivalent fractions
///
/// @eqc{\beta_{i}=-\dfrac{x_{i}z_{\mathrm{e},i}}{\sum_{j}x_{j}z_{\mathrm{e},j}},}
///
/// where @eq{x_{i}} and @eq{x_{j}} are species mole fractions and
/// @eq{z_{\mathrm{e},i}}} and @eq{z_{\mathrm{e},j}}} are exchanger
/// equivalents (or cation charges) in ion exchange species, or
/// the molar fractions
///
/// @eqc{\beta^M_{i}=-\dfrac{x_{i}}{\sum_{j}x_{j}}.}
//=====================================================================================================================

} // namespace Reaktoro
