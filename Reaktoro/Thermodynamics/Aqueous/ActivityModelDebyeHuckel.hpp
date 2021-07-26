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

/// The parameters in the Debye--Hückel activity model for aqueous solutions.
/// @ingroup ActivityModels
struct ActivityModelDebyeHuckelParams
{
    /// The default value of the *å* parameter for ionic species.
    real aiondefault = 0.0;

    /// The default value of the *b* parameter for ionic species.
    real biondefault = 0.0;

    /// The default value of the *b* parameter for neutral species.
    real bneutraldefault = 0.1;

    /// The parameters *å* of specific ionic species.
    Map<String, real> aions;

    /// The parameters *b* of specific ionic species.
    Map<String, real> bions;

    /// The parameters *b* of specific neutral species.
    Map<String, real> bneutrals;

    /// Return the *å* parameter of the ionic species with given formula.
    auto aion(const ChemicalFormula& ion) const -> real;

    /// Return the *b* parameter of the ionic species with given formula.
    auto bion(const ChemicalFormula& ion) const -> real;

    /// Return the *b* parameter of the neutral species with given formula.
    auto bneutral(const ChemicalFormula& neutral) const -> real;

    /// Set the parameters *å* and *b* of the ionic species to zero.
    auto setLimitingLaw() -> void;

    /// Set the parameters *å* of the ionic species according to Kielland (1937).
    auto setKielland() -> void;

    /// Set the parameters *å* and *b* of the ionic species according to WATEQ4F.
    auto setWATEQ4F() -> void;

    /// Set the parameters *å* and *b* of the ionic species according to PHREEQC.
    auto setPHREEQC() -> void;
};

/// Return the activity model for aqueous phases based on the Debye--Hückel model.
/// @ingroup ActivityModels
auto ActivityModelDebyeHuckel() -> ActivityModelGenerator;

/// Return the activity model for aqueous phases based on the Debye--Hückel model with given custom parameters.
/// @ingroup ActivityModels
auto ActivityModelDebyeHuckel(ActivityModelDebyeHuckelParams params) -> ActivityModelGenerator;

/// Return the activity model for aqueous phases based on the Debye--Hückel limiting law model.
/// @ingroup ActivityModels
auto ActivityModelDebyeHuckelLimitingLaw() -> ActivityModelGenerator;

/// Return the activity model for aqueous phases based on the Debye--Hückel model with Kielland (1937) parameters.
/// @ingroup ActivityModels
auto ActivityModelDebyeHuckelKielland() -> ActivityModelGenerator;

/// Return the activity model for aqueous phases based on the Debye--Hückel model using PHREEQC parameters.
/// @ingroup ActivityModels
auto ActivityModelDebyeHuckelPHREEQC() -> ActivityModelGenerator;

/// Return the activity model for aqueous phases based on the Debye--Hückel model using WATEQ4F parameters.
/// @ingroup ActivityModels
auto ActivityModelDebyeHuckelWATEQ4F() -> ActivityModelGenerator;


//=====================================================================================================================
/// @page ActivityModelDebyeHuckel Debye--Hückel activity model
/// The Debye--Hückel activity model for aqueous solutions.
/// An instance of this class can be used to control how activity coefficients
/// of ionic and neutral species, @eq{\gamma_i} and @eq{\gamma_n} respectively,
/// as well as the activity of solvent water,
/// @eq{a_\mathsf{H_2O(l)}}, are calculated using the Debye--Hückel activity
/// model.
///
/// The activity coefficients of \bold{ionics species} are calculated using the
/// following *modified Debye--Hückel equation*\supcite{Langmuir1997}:
///
/// \eqc{\log\gamma_{i}=-\dfrac{AZ_{i}^{2}\sqrt{I}}{1+B\mathring{a}_{i}\sqrt{I}}+b_{i}I,}
///
/// while the activity coefficients of \bold{neutral species} are calculated
/// using:
///
/// \eqc{\log\gamma_{n}=b_{n}I.}
///
/// In these equations, \eq{Z_i} is the electrical charge of the ionic species;
/// \eq{\mathring{a}_i} is the size or an effective diameter of the ionic
/// species (in \eq{\mathrm{Å}}, where \eq{1 \mathrm{Å}=10^{-10}\text{m}});
/// \eq{I} is the ionic strength of the aqueous solution (in molality),
/// calculated using:
///
/// \eqc{I=\frac{1}{2}\sum_{j}m_{j}Z_{j}^{2},}
///
/// with \eq{m_{j}} denoting the molality of the \eq{j}th ion. The constants
/// \eq{A} and \eq{B} in the Debye--Hückel model are calculated using (see
/// Anderson and Crerar (1993)\supcite{Anderson1993}, page 439, and Langmuir
/// (1997)\supcite{Langmuir1997}, page 128):
///
/// \eqc{A=1.824829238\cdot10^{6}\rho_{\mathrm{H_{2}O}}^{1/2}(\epsilon_{\mathrm{H_{2}O}}T)^{-3/2}}
///
/// and
///
/// \eqc{B=50.29158649\rho_{\mathrm{H_{2}O}}^{1/2}(\epsilon_{\mathrm{H_{2}O}}T)^{-1/2},}
///
/// with \eq{A} in \eq{\mathrm{(mol/kg)^{-1/2}}} and \eq{B} in units of
/// \eq{\mathrm{(mol/kg)^{-1/2}}/\mathrm{Å}}. In these equations, \eq{T} is
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
/// The activity of water is calculated using the following equation:
///
/// \eqc{\log
/// a_{w}=-\frac{1}{n_{w}^{\circ}}\left[\frac{m_{\Sigma}}{2.303}+\sum_{i}^{{\scriptscriptstyle
/// \mathrm{ions}}}m_{i}\log\gamma_{i}+\frac{2}{3}AI^{\frac{3}{2}}\sum_{i}^{{\scriptscriptstyle
/// \mathrm{ions}}}\sigma(\Lambda_{i})-I^{2}\sum_{i}^{{\scriptscriptstyle
/// \mathrm{ions}}}\frac{b_{i}}{Z_{i}^{2}}\right],}
///
/// which is thermodynamically consistent with the previous equations for the
/// activity coefficients of both ionic and neutral species, since it was
/// derived from the *Gibbs--Duhem equation*. In this equation,
/// \eq{n_{w}^{\circ}=55.508472} is the number of moles of water per kilogram;
/// \eq{m_{\Sigma}} is the sum of the molalities of all solutes (both ionic and
/// neutral species); and \eq{\sigma(\Lambda_{i})} is defined as:
///
/// \eqc{\sigma(\Lambda_{i})=\frac{3}{(\Lambda_{i}-1)^{3}}\left[(\Lambda_{i}-1)(\Lambda_{i}-3)+2\ln
/// \Lambda_{i}\right],}
///
/// with \eq{\Lambda_{i}} given by:
///
/// \eqc{\Lambda_{i}=1+B\mathring{a}_{i}\sqrt{I}.}
//=====================================================================================================================


//=====================================================================================================================
/// @fn auto ActivityModelDebyeHuckel() -> ActivityModelGenerator;
/// The activity model for aqueous phases based on the Debye--Hückel model.
/// @note This method is equivalent to ActivityModelDebyeHuckelPHREEQC().
//=====================================================================================================================


//=====================================================================================================================
/// @fn auto ActivityModelDebyeHuckel(const ActivityModelDebyeHuckelParams& params) -> ActivityModelGenerator;
/// The activity model for aqueous phases based on the Debye--Hückel model with given custom parameters.
//=====================================================================================================================


//=====================================================================================================================
/// @fn auto ActivityModelDebyeHuckelLimitingLaw() -> ActivityModelGenerator;
/// The activity model for aqueous phases based on the Debye--Hückel limiting law model.
/// Use this method to indicate that the activity coefficients of the ionic
/// species are calculated using the Debye--Hückel limiting law equation. In
/// this model, the Debye--Hückel parameters *å* and *b* of the ionic species
/// are zero.
//=====================================================================================================================


//=====================================================================================================================
/// @fn auto ActivityModelDebyeHuckelKielland() -> ActivityModelGenerator;
/// The activity model for aqueous phases based on the Debye--Hückel model with Kielland (1937) parameters.
/// In this model, the ion-size parameters *å* are taken from Kielland
/// (1937)\sup{\cite Kielland1937}:
///
/// | Ion                                                                                              | *å* (Ångström)
/// | -------------------------------------------------------------------------------------------------| --------------
/// | `H+`                                                                                             | 9
/// | `Li+`                                                                                            | 6
/// | `Rb+`, `Cs+`, `NH4+`, `Tl+`, `Ag+`                                                               | 2.5
/// | `K+`, `Cl-`, `Br-`, `I-`, `CN-`, `NO2-`, `NO3-`                                                  | 3
/// | `OH-`, `F-`, `NCS-`, `NCO-`, `HS-`, `ClO3-`, `ClO4-`, `BrO3-`, `IO4-`, `MnO4-`                   | 3.5
/// | `Na+`, `CdCl+`, `ClO2-`, `IO3-`, `HCO3-`, `H2PO4-`, `HSO3-`, `H2AsO4-`, `Co(NH3)4(NO2)2+`        | 4-4.5
/// | `Hg2+2`, `SO4-2`, `S2O3-2`, `S2O6-2`, `S2O8-2`, `SeO4-2`, `CrO4-2`, `HPO4-2`                     | 4
/// | `Pb+2`, `CO3-2`, `SO3-2`, `MoO4-2`, `Co(NH3)5Cl+2`, `Fe(CN)5NO-2`                                | 4.5
/// | `Sr+2`, `Ba+2`, `Ra+2`, `Cd+2`, `Hg+2`, `S-2`, `S2O4-2`, `WO4-2`                                 | 5
/// | `Ca+2`, `Cu+2`, `Zn+2`, `Sn+2`, `Mn+2`, `Fe+2`, `Ni+2`, `Co+2`                                   | 6
/// | `Mg+2`, `Be+2`                                                                                   | 8
/// | `PO4-3`, `Fe(CN)6-3`, `Cr(NH3)6+3`, `Co(NH3)6+3`, `Co(NH3)5H2O+3`                                | 4
/// | `Al+3`, `Fe+3`, `Cr+3`, `Sc+3`, `Y+3`, `La+3`, `In+3`, `Ce+3`, `Pr+3`, `Nd+3`, `Sm+3`            | 9
/// | `Fe(CN)6-4`                                                                                      | 5
/// | `Co(S2O3)(CN)5-4`                                                                                | 6
/// | `Th+4`, `Zn+4`, `Ce+4`, `Sn+4`                                                                   | 11
/// | `Co(SO3)2(CN)4-5`                                                                                | 9
//=====================================================================================================================


//=====================================================================================================================
/// @fn auto ActivityModelDebyeHuckelPHREEQC() -> ActivityModelGenerator;
/// The activity model for aqueous phases based on the Debye--Hückel model using PHREEQC parameters.
/// This method sets the ion-size parameters *å* and the parameter *b* of
/// the ionic species according to those used in PHREEQC v3. Their values
/// were taken from the database file `phreeqc.dat` and are listed below:
///
/// | Ion            | *å* (Å) | *b*   | Ion            | *å* (Å) | *b*
/// | -              | -       | -     | -              | -       | -
/// | `Al(OH)2+`     | 5.4     | 0     | `Al(OH)4-`     | 4.5     | 0
/// | `Al(SO4)2-`    | 4.5     | 0     | `Al+++`        | 9       | 0
/// | `AlF++`        | 5.4     | 0     | `AlF2+`        | 5.4     | 0
/// | `AlF4-`        | 4.5     | 0     | `AlOH++`       | 5.4     | 0
/// | `AlSO4+`       | 4.5     | 0     | `Ba++`         | 4       | 0.153
/// | `BaOH+`        | 5       | 0     | `Br-`          | 3       | 0
/// | `CO3--`        | 5.4     | 0     | `Ca++`         | 5       | 0.165
/// | `CaH2PO4+`     | 5.4     | 0     | `CaHCO3+`      | 6       | 0
/// | `CaPO4-`       | 5.4     | 0     | `Cl-`          | 3.63    | 0.017
/// | `Cu+`          | 2.5     | 0     | `Cu++`         | 6       | 0
/// | `CuCl+`        | 4       | 0     | `CuCl2-`       | 4       | 0
/// | `CuCl3-`       | 4       | 0     | `CuCl3--`      | 5       | 0
/// | `CuCl4--`      | 5       | 0     | `CuOH+`        | 4       | 0
/// | `F-`           | 3.5     | 0     | `Fe(OH)2+`     | 5.4     | 0
/// | `Fe(OH)3-`     | 5       | 0     | `Fe(OH)4-`     | 5.4     | 0
/// | `Fe++`         | 6       | 0     | `Fe+++`        | 9       | 0
/// | `FeCl++`       | 5       | 0     | `FeCl2+`       | 5       | 0
/// | `FeF++`        | 5       | 0     | `FeF2+`        | 5       | 0
/// | `FeH2PO4+`     | 5.4     | 0     | `FeH2PO4++`    | 5.4     | 0
/// | `FeHPO4+`      | 5       | 0     | `FeOH+`        | 5       | 0
/// | `FeOH++`       | 5       | 0     | `FeSO4+`       | 5       | 0
/// | `H+`           | 9       | 0     | `H2PO4-`       | 5.4     | 0
/// | `H2SiO4--`     | 5.4     | 0     | `H3SiO4-`      | 4       | 0
/// | `HCO3-`        | 5.4     | 0     | `HPO4--`       | 5       | 0
/// | `HS-`          | 3.5     | 0     | `K+`           | 3.5     | 0.015
/// | `KHPO4-`       | 5.4     | 0     | `KSO4-`        | 5.4     | 0
/// | `Li+`          | 6       | 0     | `LiSO4-`       | 5       | 0
/// | `Mg++`         | 5.5     | 0.2   | `MgF+`         | 4.5     | 0
/// | `MgH2PO4+`     | 5.4     | 0     | `MgHCO3+`      | 4       | 0
/// | `MgOH+`        | 6.5     | 0     | `MgPO4-`       | 5.4     | 0
/// | `Mn(OH)3-`     | 5       | 0     | `Mn++`         | 6       | 0
/// | `Mn+++`        | 9       | 0     | `MnCl+`        | 5       | 0
/// | `MnCl3-`       | 5       | 0     | `MnF+`         | 5       | 0
/// | `MnHCO3+`      | 5       | 0     | `MnOH+`        | 5       | 0
/// | `NH4+`         | 2.5     | 0     | `NO2-`         | 3       | 0
/// | `NO3-`         | 3       | 0     | `Na+`          | 4.08    | 0.082
/// | `NaHPO4-`      | 5.4     | 0     | `NaSO4-`       | 5.4     | 0
/// | `OH-`          | 3.5     | 0     | `PO4---`       | 4       | 0
/// | `S--`          | 5       | 0     | `SO4--`        | 5       | -0.04
/// | `SiF6--`       | 5       | 0     | `Sr++`         | 5.26    | 0.121
/// | `SrHCO3+`      | 5.4     | 0     | `SrOH+`        | 5       | 0
/// | `Zn++`         | 5       | 0     | `ZnCl+`        | 4       | 0
/// | `ZnCl3-`       | 4       | 0     | `ZnCl4--`      | 5       | 0
///
/// @note This method also sets the default value of *b* for neutral
/// species to 0.1, which is the default value used in PHREEQC.
///
/// **References:**
/// - Parkhurst, D. L., Appelo, C. A. J. (2013). Description of input and
///   examples for PHREEQC version 3 --- A computer program for speciation,
///   batch-reaction, one-dimensional transport, and inverse geochemical
///   calculations. In Groundwater Book 6, Modeling Techniques (p. 497).
///   U.S. Geological Survey Techniques and Methods.
//=====================================================================================================================


//=====================================================================================================================
/// @fn auto ActivityModelDebyeHuckelWATEQ4F() -> ActivityModelGenerator;
/// The activity model for aqueous phases based on the Debye--Hückel model using WATEQ4F parameters.
/// This method sets both *å* and *b* of ionic species according to the
/// ones used in WATEQ4F (Ball and Nordstrom \cite Ball1991, Truesdell and
/// Jones \cite Truesdell1974), which are listed in the following table:
///
/// | Ion          | *å*  (Å)       | *b*
/// |--------------|----------------|--------
/// | `Ca+2`       |  5.00          |  0.165
/// | `Mg+2`       |  5.50          |  0.200
/// | `Na+`        |  4.00          |  0.075
/// | `K+`         |  3.50          |  0.015
/// | `Cl-`        |  3.50          |  0.015
/// | `SO4-2`      |  5.00          | -0.040
/// | `HCO3-`      |  5.40          |  0.000
/// | `CO3-2`      |  5.40          |  0.000
/// | `Sr+2`       |  5.26          |  0.121
/// | `H+`         |  9.00          |  0.000
/// | `OH-`        |  3.50          |  0.000
/// | `SrHCO3+`    |  5.40          |  0.000
/// | `SrOH+`      |  5.00          |  0.000
/// | `Cu(S4)2-3`  | 23.00          |  0.000
/// | `CuS4S5-3`   | 25.00          |  0.000
/// | `S2-2`       |  6.50          |  0.000
/// | `S3-2`       |  8.00          |  0.000
/// | `S4-2`       | 10.00          |  0.000
/// | `S5-2`       | 12.00          |  0.000
/// | `S6-2`       | 14.00          |  0.000
/// | `Ag(S4)2-3`  | 22.00          |  0.000
/// | `AgS4S5-3`   | 24.00          |  0.000
/// | `Ag(HS)S4-2` | 15.00          |  0.000
///
/// These values for *å* and *b* are empirical. They were determined by
/// fitting the modified Debye--Hückel equation to experimental mean-salt
/// activity coefficient data.
///
/// **References:**
/// - Ball, J. W., Nordstrom, D. K. (1991). User’s Manual for WATEQ4F, with
///   revised thermodynamic data base and test cases for calculating
///   speciation of major, trace, and redox elements in natural waters.
///   U.S. Geological Survey Water-Resources Investigations Report, 91–183,
///   1–188.
/// - Truesdell, A. H., Jones, B. F. (1974). WATEQ--A computer program for
///   calculating chemical equilibrium of natural waters. U.S. Geological
///   Survey, Journal of Research, 2(2), 233–248.
//=====================================================================================================================

} // namespace Reaktoro
