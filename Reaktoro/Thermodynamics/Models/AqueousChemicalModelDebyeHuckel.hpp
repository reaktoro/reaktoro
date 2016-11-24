// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousMixture;

/// Return an equation of state for an aqueous phase based on the Debye-Huckel model.
/// @param mixture The aqueous mixture
/// @return The equation of state function for the aqueous phase
/// @see AqueousMixture, PhaseChemicalModel
auto aqueousChemicalModelDebyeHuckel(const AqueousMixture& mixture) -> PhaseChemicalModel;

/// A class used to define the parameters in the Debye--Hückel model for aqueous mixtures.
/// An instance of this class can be used to control how activity coefficients of solute
/// species, @f$\gamma_i@f$, and the activity of solvent water, @f$a_\mathsf{H_2O(l)}@f$
/// are calculated.
///
/// The activity coefficients of solute species are calculated using the following modified
/// Debye--Hückel equation (also known as B-dot equation or WATEQ Debye--Hückel equation):
/// @f[
/// \log\gamma_{i}=-\dfrac{AZ_{i}^{2}\sqrt{I}}{1+B\mathring{a}_{i}\sqrt{I}}+b_{i}I,
/// @f]
/// where @f$ Z_i @f$ is the electrical charge of the species, which is zero for neutral species;
/// @f$ I @f$ is the ionic strength of the aqueous solution (in molality); @f$ A @f$ and @f$ B @f$
/// are Debye--Hückel parameters ca
/// @f$ \mathring{a}_{i} @f$ is an ion-size parameter, which is not needed for neutral species;
///
/// species,
/// set the parameters for the modified Debye--Hückel
/// equation:
///
/// Debye--Hückel model
/// extended Debye--Hückel equation:
/// @f[
/// \log\gamma_{i}=-\dfrac{Az_{i}^{2}\sqrt{I}}{1+B\mathring{a}_{i}\sqrt{I}}
/// @f]
/// is used for those ions with known ion-size parameter @f$ \mathring{a}_{i} @f$, while the
/// Davies equation:
/// @f[
/// \log\gamma_{i}=-Az_{i}^{2}\left(\dfrac{\sqrt{I}}{1+\sqrt{I}}-0.3I\right)
/// @f]
/// is used for those ions with unknown @f$ \mathring{a}_{i} @f$.
/// An instance of this class can be used to set the parameters for the modified Debye--Hückel
/// equation:
/// @f[
/// \log\gamma_{i}=-\dfrac{Az_{i}^{2}\sqrt{I}}{1+B\mathring{a}_{i}\sqrt{I}}+b_{i}I.
/// @f]
/// Use member methods @ref a and @ref b to set the effective ion-size parameter,
/// @f$ \mathring{a}_{i} @f$ and the parameter @f$ b_{i} @f$ in the equation above.
///
/// ~~~
/// DebyeHuckel debyehuckel;
/// debyehuckel.a("Ca++") = 5.0;
/// debyehuckel.b("CO2(aq)") = 5.0;
///
/// Use the method @ref setLimitingLaw to set @f$\mathring{a}_{i}@f$ and @f$ b_{i} @f$ to zero for
/// all species, which reduces the Debye--Hückel equation to its *limiting law* form:
/// @f[
/// \log\gamma_{i}=-Az_{i}^{2}\sqrt{I}.
/// @f]
///
/// Use the method @ref setTruesdellJones to adjust the ion-size parameters @f$\mathring{a}_{i}@f$
/// and the add-on according to the table below:
///
/// | Ion        | a (Ångström) | b
/// |------------|--------------|--------
/// | Ca++       | 5.0          | 0.165
/// | Mg++       | 5.5          | 0.20
/// | Na+        | 4.0          | 0.075
/// | K+         | 3.5          | 0.015
/// | Cl-        | 3.5          | 0.015
/// | SO4--      | 5.0          | -0.04
/// | HCO3-      | 5.4          | 0.0
/// | CO3--      | 5.4          | 0.0
/// | H2CO3(aq)  | 0.0          | 0.0
/// | Sr++       | 5.26         | 0.121
/// | H+         | 9.0          | 0.0
/// | OH-        | 3.5          | 0.0
/// | SrHCO3+    | 5.4          | 0.0
/// | SrOH+      | 5.0          | 0.0
/// | SrCO3(aq)  | 0.0          | 0.0
/// | Cu(S4)2--- | 23.0         | 0.0
/// | CuS4S5---  | 25.0         | 0.0
/// | S2--       | 6.5          | 0.0
/// | S3--       | 8.0          | 0.0
/// | S4--       | 10.0         | 0.0
/// | S5--       | 12.0         | 0.0
/// | S6--       | 14.0         | 0.0
/// | Ag(S4)2--- | 22.0         | 0.0
/// | AgS4S5---  | 24.0         | 0.0
/// | Ag(HS)S4-- | 15.0         | 0.0
///
///
/// | Ion                                                                        | a (Ångström)
/// | ---------------------------------------------------------------------------|--------------
/// | H+                                                                         | 9
/// | Li+                                                                        | 6
/// | Rb+, Cs+, NH4+, Tl+, Ag+                                                   | 2.5
/// | K+, Cl-, Br-, I-, CN-, NO2-, NO3-                                          | 3
/// | OH-, F-, NCS-, NCO-, HS-, ClO3-, ClO4-, BrO3-, IO4-, MnO4-                 | 3.5
/// | Na+, CdCl+, ClO2-, IO3-, HCO3-, H2PO4-, HSO3-, H2AsO4-, Co(NH3)4(NO2)2+    | 4-4.5
/// | Hg2++, SO4--, S2O3--, S2O6--, S2O8--, SeO4--, CrO4--, HPO4--               | 4
/// | Pb++, CO3--, SO3--, MoO4--, Co(NH3)5Cl++, Fe(CN)5NO--                      | 4.5
/// | Sr++, Ba++, Ra++, Cd++, Hg++, S--, S2O4--, WO4--                           | 5
/// | Ca++, Cu++, Zn++, Sn++, Mn++, Fe++, Ni++, Co++                             | 6
/// | Mg++, Be++                                                                 | 8
/// | PO4---, Fe(CN)6---, Cr(NH3)6+++, Co(NH3)6+++, Co(NH3)5H2O+++               | 4
/// | Co(ethylenediamine)3+++                                                    | 6
/// | Al+++, Fe+++, Cr+++, Sc+++, Y+++, La+++, In+++, Ce+++, Pr+++, Nd+++, Sm+++ | 9
/// | Fe(CN)6----                                                                | 5
/// | Co(S2O3)(CN)5----                                                          | 6
/// | Th++++, Zn++++, Ce++++, Sn++++                                             | 11
/// | Co(SO3)2(CN)4-----                                                         | 9
///
///
/// | Species        | a     | b     | Species        | a     | b     | Species        | a     | b     | Species        | a     | b
/// | -              | -     | -     | -              | -     | -     | -              | -     | -     | -              | -     | -
/// | `Al(OH)2+`     | 5.4   | 0     | `Al(OH)4-`     | 4.5   | 0     | `Al(SO4)2-`    | 4.5   | 0     | `Al+++`        | 9     | 0
/// | `AlF++`        | 5.4   | 0     | `AlF2+`        | 5.4   | 0     | `AlF4-`        | 4.5   | 0     | `AlOH++`       | 5.4   | 0
/// | `AlSO4+`       | 4.5   | 0     | `Ba++`         | 4     | 0.153 | `BaOH+`        | 5     | 0     | `Br-`          | 3     | 0
/// | `CO3--`        | 5.4   | 0     | `Ca++`         | 5     | 0.165 | `CaH2PO4+`     | 5.4   | 0     | `CaHCO3+`      | 6     | 0
/// | `CaPO4-`       | 5.4   | 0     | `Cl-`          | 3.63  | 0.017 | `Cu+`          | 2.5   | 0     | `Cu++`         | 6     | 0
/// | `CuCl+`        | 4     | 0     | `CuCl2-`       | 4     | 0     | `CuCl3-`       | 4     | 0     | `CuCl3--`      | 5     | 0
/// | `CuCl4--`      | 5     | 0     | `CuOH+`        | 4     | 0     | `F-`           | 3.5   | 0     | `Fe(OH)2+`     | 5.4   | 0
/// | `Fe(OH)3-`     | 5     | 0     | `Fe(OH)4-`     | 5.4   | 0     | `Fe++`         | 6     | 0     | `Fe+++`        | 9     | 0
/// | `FeCl++`       | 5     | 0     | `FeCl2+`       | 5     | 0     | `FeF++`        | 5     | 0     | `FeF2+`        | 5     | 0
/// | `FeH2PO4+`     | 5.4   | 0     | `FeH2PO4++`    | 5.4   | 0     | `FeHPO4+`      | 5     | 0     | `FeOH+`        | 5     | 0
/// | `FeOH++`       | 5     | 0     | `FeSO4+`       | 5     | 0     | `H+`           | 9     | 0     | `H2PO4-`       | 5.4   | 0
/// | `H2SiO4--`     | 5.4   | 0     | `H3SiO4-`      | 4     | 0     | `HCO3-`        | 5.4   | 0     | `HPO4--`       | 5     | 0
/// | `HS-`          | 3.5   | 0     | `K+`           | 3.5   | 0.015 | `KHPO4-`       | 5.4   | 0     | `KSO4-`        | 5.4   | 0
/// | `Li+`          | 6     | 0     | `LiSO4-`       | 5     | 0     | `Mg++`         | 5.5   | 0.2   | `MgF+`         | 4.5   | 0
/// | `MgH2PO4+`     | 5.4   | 0     | `MgHCO3+`      | 4     | 0     | `MgOH+`        | 6.5   | 0     | `MgPO4-`       | 5.4   | 0
/// | `Mn(OH)3-`     | 5     | 0     | `Mn++`         | 6     | 0     | `Mn+++`        | 9     | 0     | `MnCl+`        | 5     | 0
/// | `MnCl3-`       | 5     | 0     | `MnF+`         | 5     | 0     | `MnHCO3+`      | 5     | 0     | `MnOH+`        | 5     | 0
/// | `NH4+`         | 2.5   | 0     | `NO2-`         | 3     | 0     | `NO3-`         | 3     | 0     | `Na+`          | 4.08  | 0.082
/// | `NaHPO4-`      | 5.4   | 0     | `NaSO4-`       | 5.4   | 0     | `OH-`          | 3.5   | 0     | `PO4---`       | 4     | 0
/// | `S--`          | 5     | 0     | `SO4--`        | 5     | -0.04 | `SiF6--`       | 5     | 0     | `Sr++`         | 5.26  | 0.121
/// | `SrHCO3+`      | 5.4   | 0     | `SrOH+`        | 5     | 0     | `Zn++`         | 5     | 0     | `ZnCl+`        | 4     | 0
/// | `ZnCl3-`       | 4     | 0     | `ZnCl4--`      | 5     | 0     |       |
///
/// **References:**
/// ----------------------------------------------------------------------------------------------
/// - Ball, J. W., Nordstrom, D. K. (1991). User’s Manual for WATEQ4F, with revised thermodynamic
///   data base and test cases for calculating speciation of major, trace, and redox elements in
///   natural waters. U.S. Geological Survey Water-Resources Investigations Report, 91–183, 1–188.
/// - Kielland, J. (1937). Individual Activity Coefficients of Ions in Aqueous Solutions.
///   Journal of the American Chemical Society, 59(9), 1675–1678.
/// - Parkhurst, D. L., Appelo, C. A. J. (2013). Description of input and examples for PHREEQC
///   version 3 — A computer program for speciation, batch-reaction, one-dimensional transport, and
///   inverse geochemical calculations. In Groundwater Book 6, Modeling Techniques (p. 497).
///   U.S. Geological Survey Techniques and Methods.
/// ~~~
class DebyeHuckel
{
public:
	/// Construct a default DebyeHuckel instance.
	DebyeHuckel();

	/// Return a reference to the `å` parameter of a given aqueous species.
	/// @param species The name of the species
	auto a(std::string species) -> double&;

	/// Return the value of the `å` parameter of a given aqueous species.
	/// @param species The name of the species
	auto a(std::string species) const -> double;

	/// Set the `å` parameter of all species to a common value.
	/// @param value The common value for the `å` parameter of all species (in units of angstrom)
	auto a(double value) -> void;

	/// Return a reference to the `b` parameter of a given aqueous species.
	/// @param species The name of the species
	auto b(std::string species) -> double&;

	/// Return the value of the `b` parameter of a given aqueous species.
	/// @param species The name of the species
	auto b(std::string species) const -> double;

	/// Set the `b` parameter of all species to a common value.
	/// @param value The common value for the `å` parameter of all species (in units of angstrom)
	auto b(double value) -> void;


	auto setLimitingLaw() -> void;


	/// Set the Debye--Hückel parameters `å` and `b` of the species according to Kielland (1937).
	/// This method sets all `b` parameters of ions to zero and the ion-size parameters `å`
	/// according to Kielland (1937). Thus, the extended Debye--Hückel equation:
	/// @f[
	/// \log\gamma_{i}=-\dfrac{Az_{i}^{2}\sqrt{I}}{1+B\mathring{a}_{i}\sqrt{I}}
	/// @f]
	/// is used for those ions with known ion-size parameter @f$ \mathring{a}_{i} @f$, while the
	/// Davies equation:
	/// @f[
	/// \log\gamma_{i}=-Az_{i}^{2}\left(\dfrac{\sqrt{I}}{1+\sqrt{I}}-0.3I\right)
	/// @f]
	/// is used for those ions with unknown @f$ \mathring{a}_{i} @f$.
	/// **References:**
	/// - Kielland, J. (1937). Individual Activity Coefficients of Ions in Aqueous Solutions.
	///   Journal of the American Chemical Society, 59(9), 1675–1678.
	auto setKielland1937() -> void;

	/// Set the Debye--Hückel parameters `å` and `b` of the species according to WATEQ4F v3.
	auto setWATEQ4F() -> void;

	/// Set the Debye--Hückel parameters `å` and `b` of the species according to PHREEQC v3.
	auto setPHREEQC() -> void;

private:

	std::map<std::string, double> m_a = a_phreeqc;
	std::map<std::string, double> m_b = b_phreeqc;

	/// The Debye--Hückel parameter `å` used in WATEQ4F (Truesdell and Jones, 1974)
	const std::map<std::string, double> a_wateq4f = {{"Ca++", 5.0}, {"Mg++", 5.5}, {"Na+", 4.0}, {"K+", 3.5}, {"Cl-", 3.5}, {"SO4--", 5.0}, {"HCO3-", 5.4}, {"CO3--", 5.4}, {"Sr++", 5.26}, {"H+", 9.0}, {"OH-", 3.5}, {"SrHCO3+", 5.4}, {"SrOH+", 5.0}, {"Cu(S4)2---", 23.0}, {"CuS4S5---", 25.0}, {"S2--", 6.5}, {"S3--", 8.0}, {"S4--", 10.0}, {"S5--", 12.0}, {"S6--", 14.0}, {"Ag(S4)2---", 22.0}, {"AgS4S5---", 24.0}, {"Ag(HS)S4--", 15.0}};

	/// The parameter `b` of the aqueous species.
	/// The Debye--Hückel parameter `å` used in WATEQ4F (Truesdell and Jones, 1974)
	const std::map<std::string, double> b_wateq4f = {{"Ca++", 0.165}, {"Mg++", 0.20}, {"Na+", 0.075}, {"K+", 0.015}, {"Cl-", 0.015}, {"SO4--", -0.04}, {"HCO3-", 0.0}, {"CO3--", 0.0}, {"H2CO3(aq)", 0.0}, {"Sr++", 0.121}};

	/// The Debye--Hückel parameter `å` used in PHREEQC v3 (Parkhurst and Appelo, 2013)
	const std::map<std::string, double> a_phreeqc = {{"Al(OH)2+", 5.4}, {"Al(OH)4-", 4.5}, {"Al(SO4)2-", 4.5}, {"Al+++", 9}, {"AlF++", 5.4}, {"AlF2+", 5.4}, {"AlF4-", 4.5}, {"AlOH++", 5.4}, {"AlSO4+", 4.5}, {"Ba++", 4}, {"BaOH+", 5}, {"Br-", 3}, {"CO3--", 5.4}, {"Ca++", 5}, {"CaH2PO4+", 5.4}, {"CaHCO3+", 6}, {"CaPO4-", 5.4}, {"Cl-", 3.63}, {"Cu+", 2.5}, {"Cu++", 6}, {"CuCl+", 4}, {"CuCl2-", 4}, {"CuCl3-", 4}, {"CuCl3--", 5}, {"CuCl4--", 5}, {"CuOH+", 4}, {"F-", 3.5}, {"Fe(OH)2+", 5.4}, {"Fe(OH)3-", 5}, {"Fe(OH)4-", 5.4}, {"Fe++", 6}, {"Fe+++", 9}, {"FeCl++", 5}, {"FeCl2+", 5}, {"FeF++", 5}, {"FeF2+", 5}, {"FeH2PO4+", 5.4}, {"FeH2PO4++", 5.4}, {"FeHPO4+", 5}, {"FeOH+", 5}, {"FeOH++", 5}, {"FeSO4+", 5}, {"H+", 9}, {"H2PO4-", 5.4}, {"H2SiO4--", 5.4}, {"H3SiO4-", 4}, {"HCO3-", 5.4}, {"HPO4--", 5}, {"HS-", 3.5}, {"K+", 3.5}, {"KHPO4-", 5.4}, {"KSO4-", 5.4}, {"Li+", 6}, {"LiSO4-", 5}, {"Mg++", 5.5}, {"MgF+", 4.5}, {"MgH2PO4+", 5.4}, {"MgHCO3+", 4}, {"MgOH+", 6.5}, {"MgPO4-", 5.4}, {"Mn(OH)3-", 5}, {"Mn++", 6}, {"Mn+++", 9}, {"MnCl+", 5}, {"MnCl3-", 5}, {"MnF+", 5}, {"MnHCO3+", 5}, {"MnOH+", 5}, {"NH4+", 2.5}, {"NO2-", 3}, {"NO3-", 3}, {"Na+", 4.08}, {"NaHPO4-", 5.4}, {"NaSO4-", 5.4}, {"OH-", 3.5}, {"PO4---", 4}, {"S--", 5}, {"SO4--", 5}, {"SiF6--", 5}, {"Sr++", 5.26}, {"SrHCO3+", 5.4}, {"SrOH+", 5}, {"Zn++", 5}, {"ZnCl+", 4}, {"ZnCl3-", 4}, {"ZnCl4--", 5}};

	/// The Debye--Hückel parameter `b` used in PHREEQC v3 (Parkhurst and Appelo, 2013)
	const std::map<std::string, double> b_phreeqc = {{"Ba++", 0.153}, {"Ca++", 0.165}, {"Cl-", 0.017}, {"K+", 0.015}, {"Mg++", 0.2}, {"Na+", 0.082}, {"SO4--", -0.04}, {"Sr++", 0.121}};
};

} // namespace Reaktoro
