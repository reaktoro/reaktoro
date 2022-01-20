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
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>

namespace Reaktoro {

/// The parameters for the surface complexation activity model.
/// @ingroup ActivityModels
struct ActivityModelSurfaceComplexationParams
{
    /// The complexation surface.
    ComplexationSurface surface;
};

/// The parameters in the Diffuse Double Layer (DDL) model for the diffusive double layer solution between complexation
/// surface and aqueous solution.
/// @ingroup ActivityModels
struct ActivityModelSurfaceComplexationDDLParams
{
    /// The default value of the enrichment factor in the Donnan pore space.
    real enrichment_factor_dll = 1.0;

    /// The diffusive double layer thickness.
    real thickness = 1e-8;

    /// The Debye length of the diffusive double layer.
    real debye_lengths = 1.0;

    /// The amount of water contained in diffusive layer can be limited.
    real limit = 0.8; // the fraction of the total water (pore space plus diffuse double layer water) that can be in the diffuse double layer

    /// The viscosity.
    real viscosity = 1.0;

    /// The complexation surface.
    ComplexationSurface surface;
};

/// Return the activity model for surface complexation assuming not encountering electrostatic effects.
/// @ingroup ActivityModels
auto ActivityModelSurfaceComplexationNoDDL(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator;

/// Return the activity model for surface complexation assuming not encountering electrostatic effects.
/// @ingroup ActivityModels
auto ActivityModelSurfaceComplexationGainesThomas(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator;

/// Return the activity model for surface complexation encountering electrostatic effects with the Diffuse Double Layer
/// (DDL) and assuming default parameters for it.
/// @ingroup ActivityModels
auto ActivityModelSurfaceComplexationDDL() -> ActivityModelGenerator;

/// Return the activity model for surface complexation encountering electrostatic effects with the Diffuse Double Layer
/// (DDL) and providing parameters of it.
/// @ingroup ActivityModels
auto ActivityModelSurfaceComplexationDDL(ActivityModelSurfaceComplexationDDLParams params) -> ActivityModelGenerator;

//=====================================================================================================================
/// @page PageActivityModelSurfaceComplexationDDL Diffuse Double Layer (DDL) activity model
///
/// The Diffuse Double Layer (DDL) model (from Dzombak & Morel, 1990), the Gouy-Chapman Diffuse Layer
/// Apparent equilibrium constant Kapp (the value of which is dependent on the surface charge) is calculated with
/// a *coulombic correction* used by Dzombak and Morel (1990), i.e.,
/// K_{app} = K_{int} \exp{\Delta Z F \Psi /RT},
/// where \Delta Z is the change in the charge of the surface species due to
/// sorption reactions, \Psi is the surface potential in volts, F is the Faraday
/// constant (96,485 C/mol), R is the molar gas constant (8.314 j/mol.K)
/// and T is the absolute temperature in K.
///
//=====================================================================================================================


} // namespace Reaktoro
