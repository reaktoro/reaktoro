//// Reaktoro is a unified framework for modeling chemically reactive systems.
////
//// Copyright (C) 2014-2015 Allan Leal
////
//// This program is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//// GNU General Public License for more details.
////
//// You should have received a copy of the GNU General Public License
//// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
//#include "AqueousChemicalModelHKF.hpp"
//
//// C++ includes
//#include <map>
//#include <string>
//#include <vector>
//
//// Reaktoro includes
//#include <Reaktoro/Common/ConvertUtils.hpp>
//#include <Reaktoro/Common/Index.hpp>
//#include <Reaktoro/Common/NamingUtils.hpp>
//#include <Reaktoro/Math/BilinearInterpolator.hpp>
//#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
//#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
//#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
//
//namespace Reaktoro {
//
///// A class used to define the parameters for the Debye--Hückel model for aqueous mixtures.
///// An instance of this class can be used to set the parameters for the modified Debye--Hückel
///// equation:
///// @f[
///// \log\gamma_{i}=-\dfrac{Az_{i}^{2}\sqrt{I}}{1+B\mathring{a}_{i}\sqrt{I}}+b_{i}I.
///// @f]
///// Use member methods @ref a and @ref b to set the effective ion-size parameter,
///// @f$ \mathring{a}_{i} @f$ and the parameter @f$ b_{i} @f$ in the equation above.
/////
///// ~~~
///// DebyeHuckel debyehuckel;
///// debyehuckel.a("Ca++") = 5.0;
///// debyehuckel.b("CO2(aq)") = 5.0;
/////
///// Use the method @ref setLimitingLaw to set @f$\mathring{a}_{i}@f$ and @f$ b_{i} @f$ to zero for
///// all species, which reduces the Debye--Hückel equation to its *limiting law* form:
///// @f[
///// \log\gamma_{i}=-Az_{i}^{2}\sqrt{I}.
///// @f]
/////
///// Use the method @ref setTruesdellJones to adjust the ion-size parameters @f$\mathring{a}_{i}@f$
///// and the add-on according to the table below:
/////
///// | Ion        | a (Ångström) | b
///// |------------|--------------|--------
///// | Ca++       | 5.0          | 0.165
///// | Mg++       | 5.5          | 0.20
///// | Na+        | 4.0          | 0.075
///// | K+         | 3.5          | 0.015
///// | Cl-        | 3.5          | 0.015
///// | SO4--      | 5.0          | -0.04
///// | HCO3-      | 5.4          | 0.0
///// | CO3--      | 5.4          | 0.0
///// | H2CO3(aq)  | 0.0          | 0.0
///// | Sr++       | 5.26         | 0.121
///// | H+         | 9.0          | 0.0
///// | OH-        | 3.5          | 0.0
///// | SrHCO3+    | 5.4          | 0.0
///// | SrOH+      | 5.0          | 0.0
///// | SrCO3(aq)  | 0.0          | 0.0
///// | Cu(S4)2--- | 23.0         | 0.0
///// | CuS4S5---  | 25.0         | 0.0
///// | S2--       | 6.5          | 0.0
///// | S3--       | 8.0          | 0.0
///// | S4--       | 10.0         | 0.0
///// | S5--       | 12.0         | 0.0
///// | S6--       | 14.0         | 0.0
///// | Ag(S4)2--- | 22.0         | 0.0
///// | AgS4S5---  | 24.0         | 0.0
///// | Ag(HS)S4-- | 15.0         | 0.0
/////
/////
/////
///// | Ion                                                                        | a (Ångström)
///// | ---------------------------------------------------------------------------|--------------
///// | H+                                                                         | 9
///// | Li+                                                                        | 6
///// | Rb+, Cs+, NH4+, Tl+, Ag+                                                   | 2.5
///// | K+, Cl-, Br-, I-, CN-, NO2-, NO3-                                          | 3
///// | OH-, F-, NCS-, NCO-, HS-, ClO3-, ClO4-, BrO3-, IO4-, MnO4-                 | 3.5
///// | Na+, CdCl+, ClO2-, IO3-, HCO3-, H2PO4-, HSO3-, H2AsO4-, Co(NH3)4(NO2)2+    | 4-4.5
///// | Hg2++, SO4--, S2O3--, S2O6--, S2O8--, SeO4--, CrO4--, HPO4--               | 4
///// | Pb++, CO3--, SO3--, MoO4--, Co(NH3)5Cl++, Fe(CN)5NO--                      | 4.5
///// | Sr++, Ba++, Ra++, Cd++, Hg++, S--, S2O4--, WO4--                           | 5
///// | Ca++, Cu++, Zn++, Sn++, Mn++, Fe++, Ni++, Co++                             | 6
///// | Mg++, Be++                                                                 | 8
///// | PO4---, Fe(CN)6---, Cr(NH3)6+++, Co(NH3)6+++, Co(NH3)5H2O+++               | 4
///// | Co(ethylenediamine)3+++                                                    | 6
///// | Al+++, Fe+++, Cr+++, Sc+++, Y+++, La+++, In+++, Ce+++, Pr+++, Nd+++, Sm+++ | 9
///// | Fe(CN)6----                                                                | 5
///// | Co(S2O3)(CN)5----                                                          | 6
///// | Th++++, Zn++++, Ce++++, Sn++++                                             | 11
///// | Co(SO3)2(CN)4-----                                                         | 9
/////
/////
///// | Species             | a                   | b                   | Species             | a                   | b                   | Species             | a                   | b                   | Species             | a                   | b
///// | -                   | -                   | -                   | -                   | -                   | -                   | -                   | -                   | -                   | -                   | -                   | -                    -
///// | `Al(OH)2+`          | 5.4                 | 0                   | `Al(OH)4-`          | 4.5                 | 0                   | `Al(SO4)2-`         | 4.5                 | 0                   | `Al+++`             | 9                   | 0
///// | `AlF++`             | 5.4                 | 0                   | `AlF2+`             | 5.4                 | 0                   | `AlF4-`             | 4.5                 | 0                   | `AlOH++`            | 5.4                 | 0
///// | `AlSO4+`            | 4.5                 | 0                   | `Ba++`              | 4                   | 0.153               | `BaOH+`             | 5                   | 0                   | `Br-`               | 3                   | 0
///// | `CO3--`             | 5.4                 | 0                   | `Ca++`              | 5                   | 0.165               | `CaH2PO4+`          | 5.4                 | 0                   | `CaHCO3+`           | 6                   | 0
///// | `CaPO4-`            | 5.4                 | 0                   | `Cl-`               | 3.63                | 0.017               | `Cu+`               | 2.5                 | 0                   | `Cu++`              | 6                   | 0
///// | `CuCl+`             | 4                   | 0                   | `CuCl2-`            | 4                   | 0                   | `CuCl3-`            | 4                   | 0                   | `CuCl3--`           | 5                   | 0
///// | `CuCl4--`           | 5                   | 0                   | `CuOH+`             | 4                   | 0                   | `F-`                | 3.5                 | 0                   | `Fe(OH)2+`          | 5.4                 | 0
///// | `Fe(OH)3-`          | 5                   | 0                   | `Fe(OH)4-`          | 5.4                 | 0                   | `Fe++`              | 6                   | 0                   | `Fe+++`             | 9                   | 0
///// | `FeCl++`            | 5                   | 0                   | `FeCl2+`            | 5                   | 0                   | `FeF++`             | 5                   | 0                   | `FeF2+`             | 5                   | 0
///// | `FeH2PO4+`          | 5.4                 | 0                   | `FeH2PO4++`         | 5.4                 | 0                   | `FeHPO4+`           | 5                   | 0                   | `FeOH+`             | 5                   | 0
///// | `FeOH++`            | 5                   | 0                   | `FeSO4+`            | 5                   | 0                   | `H+`                | 9                   | 0                   | `H2PO4-`            | 5.4                 | 0
///// | `H2SiO4--`          | 5.4                 | 0                   | `H3SiO4-`           | 4                   | 0                   | `HCO3-`             | 5.4                 | 0                   | `HPO4--`            | 5                   | 0
///// | `HS-`               | 3.5                 | 0                   | `K+`                | 3.5                 | 0.015               | `KHPO4-`            | 5.4                 | 0                   | `KSO4-`             | 5.4                 | 0
///// | `Li+`               | 6                   | 0                   | `LiSO4-`            | 5                   | 0                   | `Mg++`              | 5.5                 | 0.2                 | `MgF+`              | 4.5                 | 0
///// | `MgH2PO4+`          | 5.4                 | 0                   | `MgHCO3+`           | 4                   | 0                   | `MgOH+`             | 6.5                 | 0                   | `MgPO4-`            | 5.4                 | 0
///// | `Mn(OH)3-`          | 5                   | 0                   | `Mn++`              | 6                   | 0                   | `Mn+++`             | 9                   | 0                   | `MnCl+`             | 5                   | 0
///// | `MnCl3-`            | 5                   | 0                   | `MnF+`              | 5                   | 0                   | `MnHCO3+`           | 5                   | 0                   | `MnOH+`             | 5                   | 0
///// | `NH4+`              | 2.5                 | 0                   | `NO2-`              | 3                   | 0                   | `NO3-`              | 3                   | 0                   | `Na+`               | 4.08                | 0.082
///// | `NaHPO4-`           | 5.4                 | 0                   | `NaSO4-`            | 5.4                 | 0                   | `OH-`               | 3.5                 | 0                   | `PO4---`            | 4                   | 0
///// | `S--`               | 5                   | 0                   | `SO4--`             | 5                   | -0.04               | `SiF6--`            | 5                   | 0                   | `Sr++`              | 5.26                | 0.121
///// | `SrHCO3+`           | 5.4                 | 0                   | `SrOH+`             | 5                   | 0                   | `Zn++`              | 5                   | 0                   | `ZnCl+`             | 4                   | 0
///// | `ZnCl3-`            | 4                   | 0                   | `ZnCl4--`           | 5                   | 0                   |                     |
/////
///// References:
///// ----------------------------------------------------------------------------------------------
///// - Ball, J. W., Nordstrom, D. K. (1991). User’s Manual for WATEQ4F, with revised thermodynamic
/////   data base and test cases for calculating speciation of major, trace, and redox elements in
/////   natural waters. U.S. Geological Survey Water-Resources Investigations Report, 91–183, 1–188.
///// - Kielland, J. (1937). Individual Activity Coefficients of Ions in Aqueous Solutions.
/////   Journal of the American Chemical Society, 59(9), 1675–1678.
///// ~~~
//class DebyeHuckel
//{
//public:
//	/// Construct a default DebyeHuckel instance.
//	DebyeHuckel();
//
//	/// Return a reference to the `å` parameter of a given aqueous species.
//	/// @param species The name of the species
//	auto a(std::string species) -> double&;
//
//	/// Return the value of the `å` parameter of a given aqueous species.
//	/// @param species The name of the species
//	auto a(std::string species) const -> double;
//
//	/// Set the `å` parameter of all species to a common value.
//	/// @param value The common value for the `å` parameter of all species (in units of angstrom)
//	auto a(double value) -> void;
//
//	/// Return a reference to the `b` parameter of a given aqueous species.
//	/// @param species The name of the species
//	auto b(std::string species) -> double&;
//
//	/// Return the value of the `b` parameter of a given aqueous species.
//	/// @param species The name of the species
//	auto b(std::string species) const -> double;
//
//	/// Set the `b` parameter of all species to a common value.
//	/// @param value The common value for the `å` parameter of all species (in units of angstrom)
//	auto b(double value) -> void;
//
//	auto setLimitingLaw() -> void;
//
//	auto setKielland1937() -> void;
//
//	auto setTruesdellJones1974() -> void;
//
//	auto setParkhurst1990() -> void;
//
//	auto setWATEQ() -> void;
//
//	auto setPHREEQC() -> void;
//
//private:
//	///	The effective ion-size `å` of the aqueous species (in units of angstrom)
//	std::map<std::string, double> m_a = {
//	    {"Ca++"       , 5.0},
//	    {"Mg++"       , 5.5},
//	    {"Na+"        , 4.0},
//	    {"K+"         , 3.5},
//	    {"Cl-"        , 3.5},
//	    {"SO4--"      , 5.0},
//	    {"HCO3-"      , 5.4},
//	    {"CO3--"      , 5.4},
//	    {"H2CO3(aq)"  , 0.0},
//	    {"Sr++"       , 5.26},
//	    {"H+"         , 9.0},
//	    {"OH-"        , 3.5},
//	    {"SrHCO3+"    , 5.4},
//	    {"SrOH+"      , 5.0},
//	    {"Cu(S4)2---" , 23.0},
//	    {"CuS4S5---"  , 25.0},
//	    {"S2--"       , 6.5},
//	    {"S3--"       , 8.0},
//	    {"S4--"       , 10.0},
//	    {"S5--"       , 12.0},
//	    {"S6--"       , 14.0},
//	    {"Ag(S4)2---" , 22.0},
//	    {"AgS4S5---"  , 24.0},
//	    {"Ag(HS)S4--" , 15.0}
//	};
//
//	/// The parameter `b` of the aqueous species.
//	std::map<std::string, double> m_b = {
//	    {"Ca++"       , 0.165},
//	    {"Mg++"       , 0.20},
//	    {"Na+"        , 0.075},
//	    {"K+"         , 0.015},
//	    {"Cl-"        , 0.015},
//	    {"SO4--"      , -0.04},
//	    {"HCO3-"      , 0.0},
//	    {"CO3--"      , 0.0},
//	    {"H2CO3(aq)"  , 0.0},
//	    {"Sr++"       , 0.121}
//	};
//};
//
//namespace {
//
///// The effective ion-size `å` (in angstrom) for some aqueous species for the WATEQ Debye--Hückel model.
///// References:
///// - Ball, J. W., Nordstrom, D. K. (1991). User’s Manual for WATEQ4F, with revised thermodynamic
/////   data base and test cases for calculating speciation of major, trace, and redox elements in
/////   natural waters. U.S. Geological Survey Water-Resources Investigations Report, 91–183, 1–188.
//std::map<std::string, double> amap = {
//    {"Ca++"       , 5.0},
//    {"Mg++"       , 5.5},
//    {"Na+"        , 4.0},
//    {"K+"         , 3.5},
//    {"Cl-"        , 3.5},
//    {"SO4--"      , 5.0},
//    {"HCO3-"      , 5.4},
//    {"CO3--"      , 5.4},
//    {"H2CO3(aq)"  , 0.0},
//    {"Sr++"       , 5.26},
//    {"H+"         , 9.0},
//    {"OH-"        , 3.5},
//    {"SrHCO3+"    , 5.4},
//    {"SrOH+"      , 5.0},
//    {"Cu(S4)2---" , 23.0},
//    {"CuS4S5---"  , 25.0},
//    {"S2--"       , 6.5},
//    {"S3--"       , 8.0},
//    {"S4--"       , 10.0},
//    {"S5--"       , 12.0},
//    {"S6--"       , 14.0},
//    {"Ag(S4)2---" , 22.0},
//    {"AgS4S5---"  , 24.0},
//    {"Ag(HS)S4--" , 15.0}
//};
//
///// The parameter `b` for some aqueous species for the WATEQ Debye--Hückel model.
///// References:
///// - Ball, J. W., Nordstrom, D. K. (1991). User’s Manual for WATEQ4F, with revised thermodynamic
/////   data base and test cases for calculating speciation of major, trace, and redox elements in
/////   natural waters. U.S. Geological Survey Water-Resources Investigations Report, 91–183, 1–188.
//std::map<std::string, double> b = {
//    {"Ca++"       , 0.165},
//    {"Mg++"       , 0.20},
//    {"Na+"        , 0.075},
//    {"K+"         , 0.015},
//    {"Cl-"        , 0.015},
//    {"SO4--"      , -0.04},
//    {"HCO3-"      , 0.0},
//    {"CO3--"      , 0.0},
//    {"H2CO3(aq)"  , 0.0},
//    {"Sr++"       , 0.121}
//};
//
//} // namespace
//
//auto aqueousChemicalModelDebyeHuckel(const AqueousMixture& mixture) -> PhaseChemicalModel
//{
//    // The number of species in the mixture
//    const unsigned num_species = mixture.numSpecies();
//
//    // The number of charged and neutral species in the mixture
//    const unsigned num_charged_species = mixture.numChargedSpecies();
//    const unsigned num_neutral_species = mixture.numNeutralSpecies();
//
//    // The indices of the charged and neutral species
//    const Indices& icharged_species = mixture.indicesChargedSpecies();
//    const Indices& ineutral_species = mixture.indicesNeutralSpecies();
//
//    // The index of the water species
//    const Index iwater = mixture.indexWater();
//
//    // The effective electrostatic radii of the charged species
//    std::vector<double> effective_radii;
//
//    // The electrical charges of the charged species only
//    std::vector<double> charges;
//
//    // The natural log of 10
//    const double ln10 = std::log(10);
//
//    // The molar mass of water
//    const double Mw = waterMolarMass;
//
//    // Collect the effective radii of the ions
//    for(Index idx_ion : icharged_species)
//    {
//        const AqueousSpecies& species = mixture.species(idx_ion);
//        effective_radii.push_back(effectiveIonSize(species));
//        charges.push_back(species.charge());
//    }
//
//    // Define the intermediate chemical model function of the aqueous mixture
//    auto model = [=](const AqueousMixtureState& state)
//    {
//        // Auxiliary references to state variables
//        const auto& T = state.T;
//        const auto& I = state.Ie;
//        const auto& x = state.x;
//        const auto& m = state.m;
//        const auto& rho = state.rho/1000; // density in g/cm3
//        const auto& epsilon = state.epsilon;
//
//        // The molar fraction of the water species and its molar derivatives
//        const auto xw = x.row(iwater);
//
//        // The ln and log10 of water molar fraction
//        const auto ln_xw = log(xw);
//        const auto log10_xw = log10(xw);
//
//        // The square root of the ionic strength
//        const auto sqrtI = sqrt(I);
//        const auto sqrt_rho = sqrt(rho);
//        const auto sqrt_T_epsilon = sqrt(T * epsilon);
//        const auto T_epsilon = T * epsilon;
//
//        // The parameters for the HKF model
//        const auto A = 1.824829238e+6 * sqrt_rho/(T_epsilon*sqrt_T_epsilon);
//        const auto B = 50.29158649 * sqrt_rho/sqrt_T_epsilon;
//
//        // The alpha parameter used in the calculation of osmotic coefficient of water
//        const auto alpha = xw/(1.0 - xw) * log10_xw;
//
//        // The osmotic coefficient of the aqueous phase
//        ChemicalScalar phi(num_species);
//
//        // The result of the equation of state
//        PhaseChemicalModelResult res(num_species);
//
//        // Set the activity coefficients of the neutral species to
//        // water molar fraction to convert it to molality scale
////        res.ln_activity_coefficients = ln_xw;
//        res.ln_activity_coefficients = 0.0;
//
//        // Loop over all charged species in the mixture
//        for(unsigned i = 0; i < num_charged_species; ++i)
//        {
//            // The index of the current charged species in the mixture
//            const Index ispecies = icharged_species[i];
//
//            // The molality of the charged species and its molar derivatives
//            const auto mi = m.row(ispecies);
//
//            // Check if the molality of the charged species is zero
//            if(mi.val == 0.0)
//                continue;
//
//            // The electrical charge of the charged species
//            const auto z = charges[i];
//            const auto z2 = z*z;
//
//            // The effective radius of the charged species
//            const auto eff_radius = effective_radii[i];
//
//            // The Debye-Huckel ion size parameter of the current ion as computed by Reed (1982) and also in TOUGHREACT
//            const auto a = (z < 0) ?
//                2.0*(eff_radius + 1.91*std::abs(z))/(std::abs(z) + 1.0) :
//                2.0*(eff_radius + 1.81*std::abs(z))/(std::abs(z) + 1.0);
//
//            // The \Lamba parameter of the HKF activity coefficient model and its molar derivatives
//            const ChemicalScalar lambda = 1.0 + a*B*sqrtI;
//
//            // The log10 activity coefficient of the charged species (in molality scale) and its molar derivatives
//            const auto log10_gi = -(A*z2*sqrtI)/lambda;
//
//            // Set the activity coefficient of the current charged species
//            res.ln_activity_coefficients[ispecies] = log10_gi * ln10;
//
//            // Check if the molar fraction of water is one
//            if(xw != 1.0)
//            {
//                // The sigma parameter of the current ion and its molar derivatives
//                const auto sigma = 3.0/pow(a*B*sqrtI, 3) * (lambda - 1.0/lambda - 2.0*log(lambda));
//
//                // The psi contribution of the current ion and its molar derivatives
//                const auto psi = A*z2*sqrtI*sigma/3.0 + alpha;
//
//                // Update the osmotic coefficient with the contribution of the current charged species
//                phi += mi * psi;
//            }
//        }
//
//        // Loop over all neutral species in the mixture
//        for(unsigned i = 0; i < num_charged_species; ++i)
//        {
//            // The index of the current neutral species in the mixture
//            const Index ispecies = ineutral_species[i];
//
//            // Set the activity coefficient of the current neutral species
//            res.ln_activity_coefficients[ispecies] = log10_gi * ln10;
//        }
//
//
//        // Set the activities of the solutes (molality scale)
//        res.ln_activities = res.ln_activity_coefficients + log(m);
//
//        // Set the activity of water (in molar fraction scale)
//        if(xw != 1.0) res.ln_activities[iwater] = ln10 * Mw * phi;
//                 else res.ln_activities[iwater] = ln_xw;
//
//        // Set the activity coefficient of water (molar fraction scale)
//        res.ln_activity_coefficients[iwater] = res.ln_activities[iwater] - ln_xw;
//
//        // Set the activity constants of aqueous species to ln(55.508472)
//        res.ln_activity_constants = std::log(55.508472);
//
//        // Set the activity constant of water to zero
//        res.ln_activity_constants[iwater] = 0.0;
//
//        return res;
//    };
//
//    // Define the chemical model function of the aqueous mixture
//    PhaseChemicalModel f = [=](double T, double P, const Vector& n)
//    {
//        // Calculate the state of the mixture
//        const AqueousMixtureState state = mixture.state(T, P, n);
//
//        return model(state);
//    };
//
//    return f;
//}
//
//} // namespace Reaktoro


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

#include "AqueousChemicalModelHKF.hpp"

// C++ includes
#include <map>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {
namespace {

/// The effective electrostatic radii of ionic species (in units of angstrom).
/// This data was taken from Table 3 of Helgeson et al. (1981).
const std::map<std::string, double> effective_radii =
{
    {"H+"  , 3.08}, {"Fe+++", 3.46},
    {"Li+" , 1.64}, {"Al+++", 3.33},
    {"Na+" , 1.91}, {"Au+++", 3.72},
    {"K+"  , 2.27}, {"La+++", 3.96},
    {"Rb+" , 2.41}, {"Gd+++", 3.79},
    {"Cs+" , 2.61}, {"In+++", 3.63},
    {"NH4+", 2.31}, {"Ca+++", 3.44},
    {"Ag+" , 2.20}, {"F-"   , 1.33},
    {"Au+" , 2.31}, {"Cl-"  , 1.81},
    {"Cu+" , 1.90}, {"Br-"  , 1.96},
    {"Mg++", 2.54}, {"I-"   , 2.20},
    {"Sr++", 3.00}, {"OH-"  , 1.40},
    {"Ca++", 2.87}, {"HS-"  , 1.84},
    {"Ba++", 3.22}, {"NO3-" , 2.81},
    {"Pb++", 3.08}, {"HCO3-", 2.10},
    {"Zn++", 2.62}, {"HSO4-", 2.37},
    {"Cu++", 2.60}, {"ClO4-", 3.59},
    {"Cd++", 2.85}, {"ReO4-", 4.23},
    {"Hg++", 2.98}, {"SO4--", 3.15},
    {"Fe++", 2.62}, {"CO3--", 2.81},
    {"Mn++", 2.68}
};

/// Calculate the effective electrostatic radius of species (in units of A).
/// @param species The aqueous species instance of the ionic species
/// @return The effective electrostatic radius of the ionic species (in units of A)
double effectiveIonicRadius(const AqueousSpecies& species)
{
    // Find the effective ionic radius of the species in `effective_radii`.
    // Note that `species` might have a different name convention than those
    // used in `effective_radii`. Thus, we need to check if the name of given
    // `species` is an alternative to a name in `effective_radii`.
    for(auto pair : effective_radii)
        if(isAlternativeChargedSpeciesName(species.name(), pair.first))
            return pair.second;

    // The electrical charge of the species
    const double Zi = species.charge();

    // Estimated effective ionci radius of the species based on TOUGHREACT approach
    if(Zi == -1) return 1.81;        // based on Cl- value
    if(Zi == -2) return 3.00;        // based on rounded average of CO3-- and SO4-- values
    if(Zi == -3) return 4.20;        // based on estimation from straight line fit with charge
    if(Zi == +1) return 2.31;        // based on NH4+ value
    if(Zi == +2) return 2.80;        // based on rounded average for +2 species in the HKF table of effective ionic radii
    if(Zi == +3) return 3.60;        // based on rounded average for +3 species in the HKF table of effective ionic radii
    if(Zi == +4) return 4.50;        // based on estimaton using HKF eq. 142
    if(Zi <  -3) return -Zi*4.2/3.0; // based on linear extrapolation
    return Zi*4.5/4.0;               // based on linear extrapolation
}

} // namespace
auto aqueousChemicalModelDebyeHuckel(const AqueousMixture& mixture) -> PhaseChemicalModel
{
    // The number of species in the mixture
    const unsigned num_species = mixture.numSpecies();

    // The number of charged species in the mixture
    const unsigned num_charged_species = mixture.numChargedSpecies();

    // The indices of the charged species
    const Indices icharged_species = mixture.indicesChargedSpecies();

    // The index of the water species
    const Index iwater = mixture.indexWater();

    // The effective electrostatic radii of the charged species
    std::vector<double> effective_radii;

    // The electrical charges of the charged species only
    std::vector<double> charges;

    // The natural log of 10
    const double ln10 = std::log(10);

    // The molar mass of water
    const double Mw = waterMolarMass;

    // Collect the effective radii of the ions
    for(Index idx_ion : icharged_species)
    {
        const AqueousSpecies& species = mixture.species(idx_ion);
        effective_radii.push_back(effectiveIonicRadius(species));
        charges.push_back(species.charge());
    }

    // Define the intermediate chemical model function of the aqueous mixture
    auto model = [=](const AqueousMixtureState& state)
    {
        // Auxiliary references to state variables
        const auto& T = state.T;
        const auto& I = state.Ie;
        const auto& x = state.x;
        const auto& m = state.m;
        const auto& rho = state.rho/1000; // density in g/cm3
        const auto& epsilon = state.epsilon;

        // The molar fraction of the water species and its molar derivatives
        const auto xw = x.row(iwater);

        // The ln and log10 of water molar fraction
        const auto ln_xw = log(xw);
        const auto log10_xw = log10(xw);

        // The square root of the ionic strength
        const auto sqrtI = sqrt(I);
        const auto sqrt_rho = sqrt(rho);
        const auto sqrt_T_epsilon = sqrt(T * epsilon);
        const auto T_epsilon = T * epsilon;

        // The parameters for the HKF model
        const auto A = 1.824829238e+6 * sqrt_rho/(T_epsilon*sqrt_T_epsilon);
        const auto B = 50.29158649 * sqrt_rho/sqrt_T_epsilon;

        // The alpha parameter used in the calculation of osmotic coefficient of water
        const auto alpha = xw/(1.0 - xw) * log10_xw;

        // The osmotic coefficient of the aqueous phase
        ChemicalScalar phi(num_species);

        // The result of the equation of state
        PhaseChemicalModelResult res(num_species);

        // Set the activity coefficients of the neutral species to
        // water molar fraction to convert it to molality scale
        res.ln_activity_coefficients = ln_xw;

        // Loop over all charged species in the mixture
        for(unsigned i = 0; i < num_charged_species; ++i)
        {
            // The index of the charged species in the mixture
            const Index ispecies = icharged_species[i];

            // The molality of the charged species and its molar derivatives
            const auto mi = m.row(ispecies);

            // Check if the molality of the charged species is zero
            if(mi.val == 0.0)
                continue;

            // The electrical charge of the charged species
            const auto z = charges[i];
            const auto z2 = z*z;

            // The effective radius of the charged species
            const auto eff_radius = effective_radii[i];

            // The Debye-Huckel ion size parameter of the current ion as computed by Reed (1982) and also in TOUGHREACT
            const auto a = (z < 0) ?
                2.0*(eff_radius + 1.91*std::abs(z))/(std::abs(z) + 1.0) :
                2.0*(eff_radius + 1.81*std::abs(z))/(std::abs(z) + 1.0);

            // The \Lamba parameter of the HKF activity coefficient model and its molar derivatives
            const ChemicalScalar lambda = 1.0 + a*B*sqrtI;

            // The log10 activity coefficient of the charged species (in molality scale) and its molar derivatives
            const auto log10_gi = -(A*z2*sqrtI)/lambda + log10_xw;

            // Set the activity coefficient of the current charged species
            res.ln_activity_coefficients[ispecies] = log10_gi * ln10;

            // Check if the molar fraction of water is one
            if(xw != 1.0)
            {
                // The sigma parameter of the current ion and its molar derivatives
                const auto sigma = 3.0/pow(a*B*sqrtI, 3) * (lambda - 1.0/lambda - 2.0*log(lambda));

                // The psi contribution of the current ion and its molar derivatives
                const auto psi = A*z2*sqrtI*sigma/3.0 + alpha;

                // Update the osmotic coefficient with the contribution of the current charged species
                phi += mi * psi;
            }
        }

        // Set the activities of the solutes (molality scale)
        res.ln_activities = res.ln_activity_coefficients + log(m);

        // Set the activity of water (in molar fraction scale)
        if(xw != 1.0) res.ln_activities[iwater] = ln10 * Mw * phi;
                 else res.ln_activities[iwater] = ln_xw;

        // Set the activity coefficient of water (molar fraction scale)
        res.ln_activity_coefficients[iwater] = res.ln_activities[iwater] - ln_xw;

        // Set the activity constants of aqueous species to ln(55.508472)
        res.ln_activity_constants = std::log(55.508472);

        // Set the activity constant of water to zero
        res.ln_activity_constants[iwater] = 0.0;

        return res;
    };

    // Define the chemical model function of the aqueous mixture
    PhaseChemicalModel f = [=](double T, double P, const Vector& n)
    {
        // Calculate the state of the mixture
        const AqueousMixtureState state = mixture.state(T, P, n);

        return model(state);
    };

    return f;
}

} // namespace Reaktoro
