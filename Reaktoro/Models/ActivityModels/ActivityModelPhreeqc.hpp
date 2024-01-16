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
#include <Reaktoro/Core/Params.hpp>

namespace Reaktoro {

// Forward declarations
class PhreeqcDatabase;

/**
 * Return an activity model generator for aqueous phases based on the model used in PHREEQC.
 *
 * The activity model is equivalent to that used in PHREEQC for aqueous
 * solutions when thermodynamic databases such as `llnl.dat`, `phreeqc.dat`, and
 * similar are used. The approach is described in Parkhurst and Appelo (2013)@supcite{Parkhurst2013}. Below
 * a summary is provided.
 *
 * ### Activity coefficients of charged species
 *
 * The activity coefficients of charged species are calculated using the
 * WATEQ-Debye-Hückel model (Truesdell and Jones, 1974)@supcite{Truesdell1974}:
 *
 * \eqc{\log_{10}=-\dfrac{A_{\mathrm{DH}}z_{i}^{2}\sqrt{I}}{1+B_{\mathrm{DH}}\mathring{a}_{i}\sqrt{I}}+b_iI}
 *
 * Here, the WATEQ-Debye-Hückel parameters \f$\mathring{a}_{i}\f$ and \f$b_i\f$
 * are supplied through the `-gamma a b` command in the PHREEQC species
 * definition. The values of \f$A{\mathrm{DH}}\f$ and \f$B_{\mathrm{DH}}\f$ are
 * computed using water's density and dielectric constant.
 *
 * In cases where the `LLNL_AQUEOUS_MODEL_PARAMETERS` data block is available in
 * the PHREEQC database, the Daveler and Wolery (1992) model takes precedence:
 *
 * \eqc{\log_{10}=-\dfrac{A_{\mathrm{LNLL}}z_{i}^{2}\sqrt{I}}{1+B_{\mathrm{LNLL}}\mathring{a}_{i}\sqrt{I}}+\dot{B}I}
 *
 * Here, \f$A_{\mathrm{LNLL}}\f$, \f$B_{\mathrm{LNLL}}\f$, and \f$\dot{B}\f$ are
 * temperature-dependent parameters obtained through interpolation, with
 * relevant data within the `LLNL_AQUEOUS_MODEL_PARAMETERS` block.
 * Alternatively, in the absence of the above models, the Davies model is
 * employed:
 *
 * \eqc{\log_{10}=-A_{\mathrm{DH}}z_{i}^{2}\left(\dfrac{\sqrt{I}}{1+\sqrt{I}}-0.3I\right)}
 *
 * ### Activity coefficients of neutral species
 *
 * The activity coefficients of neutral species are determined by:
 *
 * \eqc{\log_{10}\gamma_{i}=b_iI}
 *
 * In this equation, \f$b_i\f$ takes a value of 0.1. However, when a neutral
 * species is defined using the `-co2_llnl_gamma` option, the Drummond model is
 * utilized:
 *
 * \eqc{\ln\gamma_{i}=\left(a_{1}+a_{2}T+\dfrac{a_{3}}{T}\right)I-(a_{4}+a_{5}T)\left(\dfrac{I}{I+1}\right)}
 *
 * The necessary coefficients for the Drummond model are supplied via the
 * `-co2_coefs` input within the `LLNL_AQUEOUS_MODEL_PARAMETERS` block in the
 * PHREEQC database.
 *
 * ### Activity of water
 *
 * The activity of water in PHREEQC is computed as follows (refer to equation 18
 * in PHREEQC's v2 User Guide):
 *
 * \eqc{a_{w}=1-0.017\sum_{i,\mathrm{solutes}}\dfrac{n_{i}}{M_{w}n_{w}}}
 *
 * Here, \f$M_{w}\f$ represents the molar mass of water (in kg/mol).
 *
 * @param db The PhreeqcDatabase object containing the parameters for the activity model.
 * @return An ActivityModelGenerator object that generates the activity model when constructing the phase.
 *
 * @see PhreeqcDatabase, ActivityModelGenerator
 */
auto ActivityModelPhreeqc(PhreeqcDatabase const& db) -> ActivityModelGenerator;

} // namespace Reaktoro
