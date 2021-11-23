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

/// Return the activity model for the surface complexation.
/// @see @ref PageActivityModelSurfaceComplexationDDL
/// @ingroup ActivityModels
auto ActivityModelSurfaceComplexation() -> ActivityModelGenerator;

/// Return the Diffuse Double Layer (DDL) activity model for surface complexation.
/// @see @ref PageActivityModelIonExchangeDDL
/// @ingroup ActivityModels
auto ActivityModelSurfaceComplexationDDL() -> ActivityModelGenerator;

//=====================================================================================================================
/// @page PageActivityModelSurfaceComplexationDDL Diffuse Double Layer (DDL) activity model
///
/// The Diffuse Double Layer (DDL) model (from Dzombak & Morel, 1990)
/// Gouy-Chapman Diffuse Layer
/// Apparent equilibrium constant Kapp (the value of which is dependent on the surface charge)
///
/// A *coulombic correction* used by Dzombak and Morel (1990) is defined
/// K_{int} = K_{app} \exp{\Delta Z F \Psi/RT},
/// where \Delta Z is the change in the charge of the surface species due to
/// sorption reactions, \Psi is the surface potential in volts, F is the Faraday
/// constant (96,485 C/mol), R is the molar gas constant (8.314 j/mol.K)
/// and T is the absolute temperature in K.
///
//=====================================================================================================================


} // namespace Reaktoro
