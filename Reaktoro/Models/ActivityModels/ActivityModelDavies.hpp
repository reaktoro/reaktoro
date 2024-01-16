// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

/// The parameters in the Davies activity model for aqueous electrolyte solutions.
/// @see @ref PageActivityModelDavies
/// @ingroup Thermodynamics
struct ActivityModelDaviesParams
{
    /// The value of the *b* parameter in the Davies activity model for charged species (default is 0.3).
    real bions = 0.3;

    /// The value of the *b* parameter in the Davies activity model for neutral species (default is 0.1).
    real bneutrals = 0.1;
};

/// Return the activity model for aqueous electrolyte phases based on the Davies model.
/// @see @ref PageActivityModelDavies
/// @ingroup Thermodynamics
auto ActivityModelDavies() -> ActivityModelGenerator;

/// Return the activity model for aqueous electrolyte phases based on the Davies model with given custom parameters.
/// @see @ref PageActivityModelDavies
/// @ingroup Thermodynamics
auto ActivityModelDavies(ActivityModelDaviesParams params) -> ActivityModelGenerator;

//=====================================================================================================================
/// @page PageActivityModelDavies Davies activity model
/// The Davies activity model for aqueous electrolyte solutions.
/// An instance of this class can be used to control how activity coefficients
/// of ionic and neutral species, @eq{\gamma_i} and @eq{\gamma_n} respectively,
/// as well as the activity of solvent water,
/// @eq{a_\mathsf{H_2O(aq)}}, are calculated using the Davies activity
/// model.
///
/// The activity coefficients of \bold{ionics species} are calculated using the
/// following *Davies equation*\supcite{Langmuir1997}:
///
/// \eqc{\log\gamma_{i}=-AZ_{i}^{2}\left(\dfrac{\sqrt{I}}{1+\sqrt{I}}-b_{\mathrm{charged}}I\right),}
///
/// while the activity coefficients of \bold{neutral species} are calculated
/// using:
///
/// \eqc{\log\gamma_{n}=b_{\mathrm{neutral}}I}
///
/// In these equations, \eq{Z_i} is the electrical charge of the ionic species;
/// \eq{b_{\mathrm{charged}}} is the Davies add-on parameter common to all charged species (default value is 0.3);
/// \eq{b_{\mathrm{neutral}}} is the Davies add-on parameter common to all neutral species (default value is 0.1);
/// \eq{I} is the ionic strength of the aqueous solution (in molality),
/// calculated using:
///
/// \eqc{I=\frac{1}{2}\sum_{j}m_{j}Z_{j}^{2},}
///
/// with \eq{m_{j}} denoting the molality of the \eq{j}th ion. The constant
/// \eq{A} in the Davies model is calculated using (see Anderson and Crerar (1993)\supcite{Anderson1993},
/// page 439, and Langmuir (1997)\supcite{Langmuir1997}, page 128):
///
/// \eqc{A=1.824829238\cdot10^{6}\rho_{\mathrm{H_{2}O}}^{1/2}(\epsilon_{\mathrm{H_{2}O}}T)^{-3/2},}
///
/// with \eq{A} in \eq{\mathrm{(mol/kg)^{-1/2}}}. In these equations, \eq{T} is
/// temperature (in K); \eq{\epsilon_{\mathrm{H_{2}O}}} is the dielectric
/// constant of pure water (dimensionless), calculated using the Johnson and
/// Norton (1991) model (see @ref waterElectroPropsJohnsonNorton); and
/// \eq{\rho_{\mathrm{H_{2}O}}} is the density of pure water (in
/// \eq{\mathrm{g/cm^{3}}}), calculated using either the equation of state of
/// Haar--Gallagher--Kell (1984)\sup{\cite Haar1984} or the equation of state
/// of Wagner and Pruss (2002)\sup{\cite Wagner2002} (see @ref
/// waterThermoPropsHGK and
/// @ref waterThermoPropsWagnerPruss).
///
/// The activity of water is calculated using the following equation derived
/// from the *Gibbs-Duhem conditions* for the species chemical potentials in the
/// aqueous phase:
///
/// \eqc{\ln a_{w}=A\left[2\left(\frac{I+2\sqrt{I}}{1+\sqrt{I}}\right)-4\ln(1+\sqrt{I})-bI^{2}\right]M_{\mathrm{H_{2}O}}\ln10-\frac{1-x_{w}}{x_{w}}.}
///
/// In this equation, \eq{\mathrm{H_{2}O}=0.018015268} is the molar mass of water in kg/mol.
//=====================================================================================================================


//=====================================================================================================================
/// @fn auto ActivityModelDavies() -> ActivityModelGenerator;
/// The activity model for aqueous electrolyte phases based on the Davies model.
/// @note This method is equivalent to ActivityModelDaviesPHREEQC().
//=====================================================================================================================


//=====================================================================================================================
/// @fn auto ActivityModelDavies(const ActivityModelDaviesParams& params) -> ActivityModelGenerator;
/// The activity model for aqueous electrolyte phases based on the Davies model with given custom parameters.
//=====================================================================================================================

} // namespace Reaktoro
