// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

// C++ includes
#include <map>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousMixture;
class DebyeHuckelParams;

/// Return an equation of state for an aqueous phase based on the Debye--Hückel activity model.
/// @param mixture The aqueous mixture instance
/// @param params The parameters for the Debye--Hückel activity model.
/// @return The activity model function for the aqueous phase
/// @see AqueousMixture, DebyeHuckelParams, ActivityPropsFn
auto aqueousChemicalModelDebyeHuckel(const AqueousMixture& mixture, const DebyeHuckelParams& params)-> ActivityPropsFn;

/**
A class used to define the parameters in the Debye--Hückel activity model for aqueous mixtures.

An instance of this class can be used to control how activity coefficients of ionic and neutral
species, @eq{\gamma_i} and @eq{\gamma_n} respectively, as well as the activity of solvent water,
@eq{a_\mathsf{H_2O(l)}}, are calculated using the Debye--Hückel activity model.

The activity coefficients of \bold{ionics species} are calculated using the following *modified
Debye--Hückel equation*\supcite{Langmuir1997}:

\eqc{\log\gamma_{i}=-\dfrac{AZ_{i}^{2}\sqrt{I}}{1+B\mathring{a}_{i}\sqrt{I}}+b_{i}I,}

while the activity coefficients of \bold{neutral species} are calculated using:

\eqc{\log\gamma_{n}=b_{n}I.}

In these equations, \eq{Z_i} is the electrical charge of the ionic species; \eq{\mathring{a}_i}
is the size or an effective diameter of the ionic species (in units of \eq{\mathrm{Å}}, where
\eq{1 \mathrm{Å}=10^{-10}\text{m}}); \eq{I} is the ionic strength of the aqueous solution
(in units of molality), calculated using:

\eqc{I=\frac{1}{2}\sum_{j}m_{j}Z_{j}^{2},}

with \eq{m_{j}} denoting the molality of the \eq{j}th ion. The constants \eq{A} and \eq{B} in the
Debye--Hückel model are calculated using (see Anderson and Crerar (1993)\supcite{Anderson1993},
page 439, and Langmuir (1997)\supcite{Langmuir1997}, page 128):

\eqc{A=1.824829238\cdot10^{6}\rho_{\mathrm{H_{2}O}}^{1/2}(\epsilon_{\mathrm{H_{2}O}}T)^{-3/2}}

and

\eqc{B=50.29158649\rho_{\mathrm{H_{2}O}}^{1/2}(\epsilon_{\mathrm{H_{2}O}}T)^{-1/2},}

with \eq{A} in units of \eq{\mathrm{(mol/kg)^{-1/2}}} and \eq{B} in units of
\eq{\mathrm{(mol/kg)^{-1/2}}/\mathrm{Å}}. In these equations, \eq{T} is
temperature (in units of K); \eq{\epsilon_{\mathrm{H_{2}O}}} is the dielectric constant of
pure water (dimensionless), calculated using the Johnson and Norton (1991) model
(see @ref waterElectroStateJohnsonNorton); and \eq{\rho_{\mathrm{H_{2}O}}} is the density of
pure water (in units of \eq{\mathrm{g/cm^{3}}}), calculated using either the equation of state
of Haar--Gallagher--Kell (1984)\sup{\cite Haar1984} or the equation of state of
Wagner and Pruss (2002)\sup{\cite Wagner2002} (see @ref waterThermoStateHGK and
@ref waterThermoStateWagnerPruss).

The activity of water is calculated using the following equation:

\eqc{\log a_{w}=-\frac{1}{n_{w}^{\circ}}\left[\frac{m_{\Sigma}}{2.303}+\sum_{i}^{{\scriptscriptstyle \mathrm{ions}}}m_{i}\log\gamma_{i}+\frac{2}{3}AI^{\frac{3}{2}}\sum_{i}^{{\scriptscriptstyle \mathrm{ions}}}\sigma(\Lambda_{i})-I^{2}\sum_{i}^{{\scriptscriptstyle \mathrm{ions}}}\frac{b_{i}}{Z_{i}^{2}}\right],}

which is thermodynamically consistent with the previous equations for the activity coefficients of
both ionic and neutral species, since it was derived from the *Gibbs--Duhem equation*.
In this equation, \eq{n_{w}^{\circ}=55.508472} is the number of moles of water per kilogram;
\eq{m_{\Sigma}} is the sum of the molalities of all solutes (both ionic and neutral species);
and \eq{\sigma(\Lambda_{i})} is defined as:

\eqc{\sigma(\Lambda_{i})=\frac{3}{(\Lambda_{i}-1)^{3}}\left[(\Lambda_{i}-1)(\Lambda_{i}-3)+2\ln \Lambda_{i}\right],}

with \eq{\Lambda_{i}} given by:

\eqc{\Lambda_{i}=1+B\mathring{a}_{i}\sqrt{I}.}

**Usage Example:**

Use member methods @ref aion and @ref bion to set the effective ion-size parameter,
\eq{\mathring{a}_i}, and the parameter \eq{b_i} in the previous equations. Use method @ref bneutral
to set the parameter \eq{b_n} of a neutral species.

~~~
using namespace Reaktoro; {delete}
DebyeHuckelParams debyehuckel;
debyehuckel.aion("Ca++", 5.0);
debyehuckel.bion("Ca++", 0.165);
debyehuckel.bneutral("O2(aq)", 0.2);
~~~
*/
class DebyeHuckelParams
{
public:
	/// Construct a default DebyeHuckelParams instance.
    DebyeHuckelParams();

	/// Set the default ion-size parameter value to be used for ionic species lacking such data.
	/// @param value The default ion-size parameter value
	auto aiondefault(double value) -> void;

	/// Return the default ion-size parameter value to be used for ionic species lacking such data.
	auto aiondefault() const -> double;

	/// Set the ion-size parameter value of a given ionic species (in units of Å).
	/// @param name The name of the ionic species
	/// @param value The ion-size parameter value (in units of Å)
	auto aion(std::string name, double value) -> void;

	/// Set the ion-size parameter values of several ionic species (in units of Å).
	/// @param pairs The pairs (name, value) for the ion-size parameters of the ionic species
	auto aion(const std::map<std::string, double>& pairs) -> void;

	/// Set the ion-size parameter of all ionic species to a common value.
	/// @warning This method **overwrites** all previously assigned ion-size parameter values.
	///          It also **overwrites** the default value of the ion-size parameter that is used
	///          for those ions lacking such data.
	/// @param value The common ion-size parameter value (in units of Å)
	auto aion(double value) -> void;

	/// Return the ion-size parameter value of a given ionic species (in units of Å).
	/// @param name The name of the ionic species
	auto aion(std::string name) const -> double;

	/// Set the default Debye--Hückel parameter *b* to be used for ionic species lacking such data.
	/// @param value The default parameter *b* value for ionic species
	auto biondefault(double value) -> void;

	/// Return the default Debye--Hückel parameter *b* to be used for ionic species lacking such data.
	auto biondefault() const -> double;

	/// Set the value of the Debye--Hückel parameter *b* of a given ionic species.
	/// @param name The name of the ionic species
	/// @param value The value of the *b* parameter
	auto bion(std::string name, double value) -> void;

    /// Set the value of the Debye--Hückel parameter *b* of several ionic species (in units of Å).
    /// @param pairs The pairs (name, value) for the parameter *b* of the ionic species
    auto bion(const std::map<std::string, double>& pairs) -> void;

	/// Set the Debye--Hückel parameter *b* of all ionic species to a common value.
	/// @warning This method **overwrites** all previously assigned values for the *b* parameter
	///          of ionic species. It also **overwrites** the default *b* value that is used
	///          for those ions lacking such data.
	/// @param value The common value for the *b* parameter of the ionic species
	auto bion(double value) -> void;

	/// Return the value of the *b* parameter of a given ionic species.
	/// @param name The name of the ionic species
	auto bion(std::string name) const -> double;

	/// Set the default value of the *b* parameter to be used neutral species lacking such data.
	/// @param value The default parameter *b* value for neutral species
	auto bneutraldefault(double value) -> void;

	/// Return the default value of the *b* parameter to be used neutral species lacking such data.
	auto bneutraldefault() const -> double;

	/// Set the value of the *b* parameter of a given neutral species.
	/// @param name The name of the neutral species
	/// @param value The value of the *b* parameter
	auto bneutral(std::string name, double value) -> void;

	/// Set the value of the *b* parameter of several neutral species.
    /// @param pairs The pairs (name, value) for the parameter *b* of the neutral species
    auto bneutral(const std::map<std::string, double>& pairs) -> void;

	/// Set the value of the *b* parameter of all neutral species to a common value.
	/// @warning This method **overwrites** all previously assigned values for the *b* parameter
	///          of neutral species. It also **overwrites** the default *b* value that is used
	///          for those neutral species lacking such data.
	/// @param value The common value for the *b* parameter of neutral species
	auto bneutral(double value) -> void;

	/// Return the value of the *b* parameter of a given neutral species.
	/// @note This method returns the default *b* value for the neutral species if the given
	///       speciesof the *b* parameter for the neutral species that was
	///       previously define
	/// @param name The name of the neutral species
	auto bneutral(std::string name) const -> double;

	/// Set the Debye--Hückel limiting law model for the ionic species.
	/// Use this method to indicate that the activity coefficients of the ionic species are
	/// calculated using the Debye--Hückel limiting law equation.
	///
	/// @warning This method sets the default values of the Debye--Hückel parameters *å* and *b*
	///          of the ionic species to zero. It also **overwrites** previously set values
	///          for *å* and *b* of the ionic species.
	///
	/// @note This method is equivalent to calling methods aion and bion with argument 0.0.
	///       For example, `params.aion(0.0); params.bion(0.0);`.
	///
	/// @see aion, bion
	auto setLimitingLaw() -> void;

	/**
	Set the ion-size parameters *å* according to the values given in Kielland (1937)\sup{\cite Kielland1937}.
	This method sets the values of ion-size parameters *å* of the ionic species according to
	those given in Kielland (1937)\sup{\cite Kielland1937}, which is summarized in the
	following table:

    | Ion                                                                                              | *å* (Ångström)
    | -------------------------------------------------------------------------------------------------| --------------
    | `H+`                                                                                             | 9
    | `Li+`                                                                                            | 6
    | `Rb+`, `Cs+`, `NH4+`, `Tl+`, `Ag+`                                                               | 2.5
    | `K+`, `Cl-`, `Br-`, `I-`, `CN-`, `NO2-`, `NO3-`                                                  | 3
    | `OH-`, `F-`, `NCS-`, `NCO-`, `HS-`, `ClO3-`, `ClO4-`, `BrO3-`, `IO4-`, `MnO4-`                   | 3.5
    | `Na+`, `CdCl+`, `ClO2-`, `IO3-`, `HCO3-`, `H2PO4-`, `HSO3-`, `H2AsO4-`, `Co(NH3)4(NO2)2+`        | 4-4.5
    | `Hg2++`, `SO4--`, `S2O3--`, `S2O6--`, `S2O8--`, `SeO4--`, `CrO4--`, `HPO4--`                     | 4
    | `Pb++`, `CO3--`, `SO3--`, `MoO4--`, `Co(NH3)5Cl++`, `Fe(CN)5NO--`                                | 4.5
    | `Sr++`, `Ba++`, `Ra++`, `Cd++`, `Hg++`, `S--`, `S2O4--`, `WO4--`                                 | 5
    | `Ca++`, `Cu++`, `Zn++`, `Sn++`, `Mn++`, `Fe++`, `Ni++`, `Co++`                                   | 6
    | `Mg++`, `Be++`                                                                                   | 8
    | `PO4---`, `Fe(CN)6---`, `Cr(NH3)6+++`, `Co(NH3)6+++`, `Co(NH3)5H2O+++`                           | 4
    | `Al+++`, `Fe+++`, `Cr+++`, `Sc+++`, `Y+++`, `La+++`, `In+++`, `Ce+++`, `Pr+++`, `Nd+++`, `Sm+++` | 9
    | `Fe(CN)6----`                                                                                    | 5
    | `Co(S2O3)(CN)5----`                                                                              | 6
    | `Th++++`, `Zn++++`, `Ce++++`, `Sn++++`                                                           | 11
    | `Co(SO3)2(CN)4-----`                                                                             | 9

	@warning This method **overwrites** previously assigned values of *å* for those ionic
	         species in the table above.

	@note This method leaves unchanged the Debye--Hückel parameters *b* of the ionic species.
	*/
	auto setKielland1937() -> void;

	/// Set the Debye--Hückel parameters *å* and *b* of the ionic species according to WATEQ4F\sup{\cite Truesdell1974 \cite Ball1991}.
	/// This method sets both *å* and *b* of ionic species according to the ones used in WATEQ4F
	/// (Ball and Nordstrom \cite Ball1991, Truesdell and Jones cite Truesdell1974), which is listed
	/// in the following table:
	///
	/// | Ion          | *å*  (Å)       | *b*
    /// |--------------|----------------|--------
    /// | `Ca++`       | 5.0            | 0.165
    /// | `Mg++`       | 5.5            | 0.20
    /// | `Na+`        | 4.0            | 0.075
    /// | `K+`         | 3.5            | 0.015
    /// | `Cl-`        | 3.5            | 0.015
    /// | `SO4--`      | 5.0            | -0.04
    /// | `HCO3-`      | 5.4            | 0.0
    /// | `CO3--`      | 5.4            | 0.0
    /// | `Sr++`       | 5.26           | 0.121
    /// | `H+`         | 9.0            | 0.0
    /// | `OH-`        | 3.5            | 0.0
    /// | `SrHCO3+`    | 5.4            | 0.0
    /// | `SrOH+`      | 5.0            | 0.0
    /// | `Cu(S4)2---` | 23.0           | 0.0
    /// | `CuS4S5---`  | 25.0           | 0.0
    /// | `S2--`       | 6.5            | 0.0
    /// | `S3--`       | 8.0            | 0.0
    /// | `S4--`       | 10.0           | 0.0
    /// | `S5--`       | 12.0           | 0.0
    /// | `S6--`       | 14.0           | 0.0
    /// | `Ag(S4)2---` | 22.0           | 0.0
    /// | `AgS4S5---`  | 24.0           | 0.0
    /// | `Ag(HS)S4--` | 15.0           | 0.0
    ///
    /// These values for *å* and *b* are empirical. They were determined by fitting the modified
    /// Debye--Hückel equation to experimental mean-salt activity coefficient data.
    ///
    /// @warning This method **overwrites** previously assigned values of *å* and *b* for those
    ///          ionic species in the table above.
    ///
	/// **References:**
	/// - Ball, J. W., Nordstrom, D. K. (1991). User’s Manual for WATEQ4F, with revised
	///   thermodynamic data base and test cases for calculating speciation of major, trace, and
	///   redox elements in natural waters. U.S. Geological Survey Water-Resources Investigations
	///   Report, 91–183, 1–188.
	/// - Truesdell, A. H., Jones, B. F. (1974). WATEQ--A computer program for calculating chemical
	///   equilibrium of natural waters. U.S. Geological Survey, Journal of Research, 2(2), 233–248.
	///
	auto setWATEQ4F() -> void;

	/// Set the Debye--Hückel parameters *å* and *b* of the species according to PHREEQC v3.
	/// This method sets the ion-size parameters *å* and the parameter *b* of the ionic species
	/// according to those used in PHREEQC v3 when the `phreeqc.dat` database file is used, which
	/// are listed in the following table:
	///
    ///	| Ion            | *å* (Å) | *b*   | Ion            | *å* (Å) | *b*
    ///	| -              | -       | -     | -              | -       | -
    ///	| `Al(OH)2+`     | 5.4     | 0     | `Al(OH)4-`     | 4.5     | 0
    ///	| `Al(SO4)2-`    | 4.5     | 0     | `Al+++`        | 9       | 0
    ///	| `AlF++`        | 5.4     | 0     | `AlF2+`        | 5.4     | 0
    ///	| `AlF4-`        | 4.5     | 0     | `AlOH++`       | 5.4     | 0
    ///	| `AlSO4+`       | 4.5     | 0     | `Ba++`         | 4       | 0.153
    ///	| `BaOH+`        | 5       | 0     | `Br-`          | 3       | 0
    ///	| `CO3--`        | 5.4     | 0     | `Ca++`         | 5       | 0.165
    ///	| `CaH2PO4+`     | 5.4     | 0     | `CaHCO3+`      | 6       | 0
    ///	| `CaPO4-`       | 5.4     | 0     | `Cl-`          | 3.63    | 0.017
    ///	| `Cu+`          | 2.5     | 0     | `Cu++`         | 6       | 0
    ///	| `CuCl+`        | 4       | 0     | `CuCl2-`       | 4       | 0
    ///	| `CuCl3-`       | 4       | 0     | `CuCl3--`      | 5       | 0
    ///	| `CuCl4--`      | 5       | 0     | `CuOH+`        | 4       | 0
    ///	| `F-`           | 3.5     | 0     | `Fe(OH)2+`     | 5.4     | 0
    ///	| `Fe(OH)3-`     | 5       | 0     | `Fe(OH)4-`     | 5.4     | 0
    ///	| `Fe++`         | 6       | 0     | `Fe+++`        | 9       | 0
    ///	| `FeCl++`       | 5       | 0     | `FeCl2+`       | 5       | 0
    ///	| `FeF++`        | 5       | 0     | `FeF2+`        | 5       | 0
    ///	| `FeH2PO4+`     | 5.4     | 0     | `FeH2PO4++`    | 5.4     | 0
    ///	| `FeHPO4+`      | 5       | 0     | `FeOH+`        | 5       | 0
    ///	| `FeOH++`       | 5       | 0     | `FeSO4+`       | 5       | 0
    ///	| `H+`           | 9       | 0     | `H2PO4-`       | 5.4     | 0
    ///	| `H2SiO4--`     | 5.4     | 0     | `H3SiO4-`      | 4       | 0
    ///	| `HCO3-`        | 5.4     | 0     | `HPO4--`       | 5       | 0
    ///	| `HS-`          | 3.5     | 0     | `K+`           | 3.5     | 0.015
    ///	| `KHPO4-`       | 5.4     | 0     | `KSO4-`        | 5.4     | 0
    ///	| `Li+`          | 6       | 0     | `LiSO4-`       | 5       | 0
    ///	| `Mg++`         | 5.5     | 0.2   | `MgF+`         | 4.5     | 0
    ///	| `MgH2PO4+`     | 5.4     | 0     | `MgHCO3+`      | 4       | 0
    ///	| `MgOH+`        | 6.5     | 0     | `MgPO4-`       | 5.4     | 0
    ///	| `Mn(OH)3-`     | 5       | 0     | `Mn++`         | 6       | 0
    ///	| `Mn+++`        | 9       | 0     | `MnCl+`        | 5       | 0
    ///	| `MnCl3-`       | 5       | 0     | `MnF+`         | 5       | 0
    ///	| `MnHCO3+`      | 5       | 0     | `MnOH+`        | 5       | 0
    ///	| `NH4+`         | 2.5     | 0     | `NO2-`         | 3       | 0
    ///	| `NO3-`         | 3       | 0     | `Na+`          | 4.08    | 0.082
    ///	| `NaHPO4-`      | 5.4     | 0     | `NaSO4-`       | 5.4     | 0
    ///	| `OH-`          | 3.5     | 0     | `PO4---`       | 4       | 0
    ///	| `S--`          | 5       | 0     | `SO4--`        | 5       | -0.04
    ///	| `SiF6--`       | 5       | 0     | `Sr++`         | 5.26    | 0.121
    ///	| `SrHCO3+`      | 5.4     | 0     | `SrOH+`        | 5       | 0
    ///	| `Zn++`         | 5       | 0     | `ZnCl+`        | 4       | 0
    ///	| `ZnCl3-`       | 4       | 0     | `ZnCl4--`      | 5       | 0
	///
	/// @note This method also sets the default value of *b* for neutral species to 0.1, which is
	///       the default value used in PHREEQC.
	///
	/// **References:**
	/// - Parkhurst, D. L., Appelo, C. A. J. (2013). Description of input and examples for PHREEQC
	///   version 3 --- A computer program for speciation, batch-reaction, one-dimensional
	///   transport, and inverse geochemical calculations. In Groundwater Book 6, Modeling
	///   Techniques (p. 497). U.S. Geological Survey Techniques and Methods.
	auto setPHREEQC() -> void;

private:
	struct Impl;

	std::shared_ptr<Impl> pimpl;
};

} // namespace Reaktoro
