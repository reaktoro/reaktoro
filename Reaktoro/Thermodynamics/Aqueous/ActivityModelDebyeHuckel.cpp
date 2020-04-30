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

#include "ActivityModelDebyeHuckel.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

using std::log;
using std::pow;
using std::sqrt;

namespace detail {

/// The Debye--Hückel parameter `å` used in PHREEQC v3 (Parkhurst and Appelo, 2013)
const Map<String, real> aions_phreeqc =
{
    { "Al(OH)2+" , 5.4  },
    { "Al(OH)4-" , 4.5  },
    { "Al(SO4)2-", 4.5  },
    { "Al+++"    , 9.0  },
    { "AlF++"    , 5.4  },
    { "AlF2+"    , 5.4  },
    { "AlF4-"    , 4.5  },
    { "AlOH++"   , 5.4  },
    { "AlSO4+"   , 4.5  },
    { "Ba++"     , 4.0  },
    { "BaOH+"    , 5.0  },
    { "Br-"      , 3.0  },
    { "CO3--"    , 5.4  },
    { "Ca++"     , 5.0  },
    { "CaH2PO4+" , 5.4  },
    { "CaHCO3+"  , 6.0  },
    { "CaPO4-"   , 5.4  },
    { "Cl-"      , 3.63 },
    { "Cu+"      , 2.5  },
    { "Cu++"     , 6.0  },
    { "CuCl+"    , 4.0  },
    { "CuCl2-"   , 4.0  },
    { "CuCl3-"   , 4.0  },
    { "CuCl3--"  , 5.0  },
    { "CuCl4--"  , 5.0  },
    { "CuOH+"    , 4.0  },
    { "F-"       , 3.5  },
    { "Fe(OH)2+" , 5.4  },
    { "Fe(OH)3-" , 5.0  },
    { "Fe(OH)4-" , 5.4  },
    { "Fe++"     , 6.0  },
    { "Fe+++"    , 9.0  },
    { "FeCl++"   , 5.0  },
    { "FeCl2+"   , 5.0  },
    { "FeF++"    , 5.0  },
    { "FeF2+"    , 5.0  },
    { "FeH2PO4+" , 5.4  },
    { "FeH2PO4++", 5.4  },
    { "FeHPO4+"  , 5.0  },
    { "FeOH+"    , 5.0  },
    { "FeOH++"   , 5.0  },
    { "FeSO4+"   , 5.0  },
    { "H+"       , 9.0  },
    { "H2PO4-"   , 5.4  },
    { "H2SiO4--" , 5.4  },
    { "H3SiO4-"  , 4.0  },
    { "HCO3-"    , 5.4  },
    { "HPO4--"   , 5.0  },
    { "HS-"      , 3.5  },
    { "K+"       , 3.5  },
    { "KHPO4-"   , 5.4  },
    { "KSO4-"    , 5.4  },
    { "Li+"      , 6.0  },
    { "LiSO4-"   , 5.0  },
    { "Mg++"     , 5.5  },
    { "MgF+"     , 4.5  },
    { "MgH2PO4+" , 5.4  },
    { "MgHCO3+"  , 4.0  },
    { "MgOH+"    , 6.5  },
    { "MgPO4-"   , 5.4  },
    { "Mn(OH)3-" , 5.0  },
    { "Mn++"     , 6.0  },
    { "Mn+++"    , 9.0  },
    { "MnCl+"    , 5.0  },
    { "MnCl3-"   , 5.0  },
    { "MnF+"     , 5.0  },
    { "MnHCO3+"  , 5.0  },
    { "MnOH+"    , 5.0  },
    { "NH4+"     , 2.5  },
    { "NO2-"     , 3.0  },
    { "NO3-"     , 3.0  },
    { "Na+"      , 4.08 },
    { "NaHPO4-"  , 5.4  },
    { "NaSO4-"   , 5.4  },
    { "OH-"      , 3.5  },
    { "PO4---"   , 4.0  },
    { "S--"      , 5.0  },
    { "SO4--"    , 5.0  },
    { "SiF6--"   , 5.0  },
    { "Sr++"     , 5.26 },
    { "SrHCO3+"  , 5.4  },
    { "SrOH+"    , 5.0  },
    { "Zn++"     , 5.0  },
    { "ZnCl+"    , 4.0  },
    { "ZnCl3-"   , 4.0  },
    { "ZnCl4--"  , 5.0  }
};

/// The Debye--Hückel parameter `b` used in PHREEQC v3 (Parkhurst and Appelo, 2013)
const Map<String, real> bions_phreeqc =
{
    {"Ba++" ,  0.153 },
    {"Ca++" ,  0.165 },
    {"Cl-"  ,  0.017 },
    {"K+"   ,  0.015 },
    {"Mg++" ,  0.200 },
    {"Na+"  ,  0.082 },
    {"SO4--", -0.040 },
    {"Sr++" ,  0.121 }
};

/// The Debye--Hückel parameter `å` used in WATEQ4F (Ball and Nordstrom 1991, Truesdell and Jones 1974)
const Map<String, real> aions_wateq4f =
{
    { "Ca++"      ,  5.00 },
    { "Mg++"      ,  5.50 },
    { "Na+"       ,  4.00 },
    { "K+"        ,  3.50 },
    { "Cl-"       ,  3.50 },
    { "SO4--"     ,  5.00 },
    { "HCO3-"     ,  5.40 },
    { "CO3--"     ,  5.40 },
    { "Sr++"      ,  5.26 },
    { "H+"        ,  9.00 },
    { "OH-"       ,  3.50 },
    { "SrHCO3+"   ,  5.40 },
    { "SrOH+"     ,  5.00 },
    { "Cu(S4)2---", 23.00 },
    { "CuS4S5---" , 25.00 },
    { "S2--"      ,  6.50 },
    { "S3--"      ,  8.00 },
    { "S4--"      , 10.00 },
    { "S5--"      , 12.00 },
    { "S6--"      , 14.00 },
    { "Ag(S4)2---", 22.00 },
    { "AgS4S5---" , 24.00 },
    { "Ag(HS)S4--", 15.00 }
};

/// The Debye--Hückel parameter `b` used in WATEQ4F (Ball and Nordstrom 1991, Truesdell and Jones 1974)
const Map<String, real> bions_wateq4f =
{
    { "Ca++"     ,  0.165 },
    { "Mg++"     ,  0.200 },
    { "Na+"      ,  0.075 },
    { "K+"       ,  0.015 },
    { "Cl-"      ,  0.015 },
    { "SO4--"    , -0.040 },
    { "HCO3-"    ,  0.000 },
    { "CO3--"    ,  0.000 },
    { "H2CO3(aq)",  0.000 },
    { "Sr++"     ,  0.121 }
};

/// The Debye--Hückel parameter `å` from Kielland (1937).
const Map<String, real> aions_kielland =
{
    { "H+"                 ,  9.0 },
    { "Li+"                ,  6.0 },
    { "Rb+"                ,  2.5 },
    { "Cs+"                ,  2.5 },
    { "NH4+"               ,  2.5 },
    { "Tl+"                ,  2.5 },
    { "Ag+"                ,  2.5 },
    { "K+"                 ,  3.0 },
    { "Cl-"                ,  3.0 },
    { "Br-"                ,  3.0 },
    { "I-"                 ,  3.0 },
    { "CN-"                ,  3.0 },
    { "NO2-"               ,  3.0 },
    { "NO3-"               ,  3.0 },
    { "OH-"                ,  3.5 },
    { "F-"                 ,  3.5 },
    { "NCS-"               ,  3.5 },
    { "NCO-"               ,  3.5 },
    { "HS-"                ,  3.5 },
    { "ClO3-"              ,  3.5 },
    { "ClO4-"              ,  3.5 },
    { "BrO3-"              ,  3.5 },
    { "IO4-"               ,  3.5 },
    { "MnO4-"              ,  3.5 },
    { "Na+"                ,  4.0 },
    { "CdCl+"              ,  4.0 },
    { "ClO2-"              ,  4.0 },
    { "IO3-"               ,  4.0 },
    { "HCO3-"              ,  4.0 },
    { "H2PO4-"             ,  4.0 },
    { "HSO3-"              ,  4.0 },
    { "H2AsO4-"            ,  4.0 },
    { "Co(NH3)4(NO2)2+"    ,  4.0 },
    { "Hg2++"              ,  4.0 },
    { "SO4--"              ,  4.0 },
    { "S2O3--"             ,  4.0 },
    { "S2O6--"             ,  4.0 },
    { "S2O8--"             ,  4.0 },
    { "SeO4--"             ,  4.0 },
    { "CrO4--"             ,  4.0 },
    { "HPO4--"             ,  4.0 },
    { "Pb++"               ,  4.5 },
    { "CO3--"              ,  4.5 },
    { "SO3--"              ,  4.5 },
    { "MoO4--"             ,  4.5 },
    { "Co(NH3)5Cl++"       ,  4.5 },
    { "Fe(CN)5NO--"        ,  4.5 },
    { "Sr++"               ,  5.0 },
    { "Ba++"               ,  5.0 },
    { "Ra++"               ,  5.0 },
    { "Cd++"               ,  5.0 },
    { "Hg++"               ,  5.0 },
    { "S--"                ,  5.0 },
    { "S2O4--"             ,  5.0 },
    { "WO4--"              ,  5.0 },
    { "Ca++"               ,  6.0 },
    { "Cu++"               ,  6.0 },
    { "Zn++"               ,  6.0 },
    { "Sn++"               ,  6.0 },
    { "Mn++"               ,  6.0 },
    { "Fe++"               ,  6.0 },
    { "Ni++"               ,  6.0 },
    { "Co++"               ,  6.0 },
    { "Mg++"               ,  8.0 },
    { "Be++"               ,  8.0 },
    { "PO4---"             ,  4.0 },
    { "Fe(CN)6---"         ,  4.0 },
    { "Cr(NH3)6+++"        ,  4.0 },
    { "Co(NH3)6+++"        ,  4.0 },
    { "Co(NH3)5H2O+++"     ,  4.0 },
    { "Al+++"              ,  9.0 },
    { "Fe+++"              ,  9.0 },
    { "Cr+++"              ,  9.0 },
    { "Sc+++"              ,  9.0 },
    { "Y+++"               ,  9.0 },
    { "La+++"              ,  9.0 },
    { "In+++"              ,  9.0 },
    { "Ce+++"              ,  9.0 },
    { "Pr+++"              ,  9.0 },
    { "Nd+++"              ,  9.0 },
    { "Sm+++"              ,  9.0 },
    { "Fe(CN)6----"        ,  5.0 },
    { "Co(S2O3)(CN)5----"  ,  6.0 },
    { "Th++++"             , 11.0 },
    { "Zn++++"             , 11.0 },
    { "Ce++++"             , 11.0 },
    { "Sn++++"             , 11.0 },
    { "Co(SO3)2(CN)4-----" ,  9.0 }
};

/// Return a parameter value whose key matches a given formula
auto get(const ChemicalFormula& formula, const Map<String, real>& params, real defaultvalue) -> real
{
    for(const auto& [key, value] : params)
        if(formula.equivalent(key))
            return value;
    return defaultvalue;
};

/// Return the ActivityPropsFn object based on the Debye-Huckel model.
auto activityPropsFnDebyeHuckel(const SpeciesList& species, ActivityModelDebyeHuckelParams params) -> ActivityPropsFn
{
    // Create the aqueous mixture
    AqueousMixture mixture(species);

    // The molar mass of water
    const auto Mw = waterMolarMass;

    // The number of moles of water per kg
    const auto nwo = 1.0/Mw;

    // The number of species in the aqueous mixture
    const auto num_species = species.size();

    // The number of charged and neutral species in the aqueous mixture
    const auto num_charged_species = mixture.charged().size();
    const auto num_neutral_species = mixture.neutral().size();

    // The indices of the charged and neutral species
    const auto icharged_species = mixture.indicesCharged();
    const auto ineutral_species = mixture.indicesNeutral();

    // The index of the water species
    const auto iwater = mixture.indexWater();

    // The electrical charges of the charged species only
    const ArrayXd charges = mixture.charges()(icharged_species);

    // The Debye-Huckel parameters a and b of the charged species
    Vec<real> aions, bions;

    // The Debye-Huckel parameter b of the neutral species
    Vec<real> bneutral;

    // Collect the Debye-Huckel parameters a and b of the charged species
    for(Index i : icharged_species)
    {
        const auto species = mixture.species(i);
        aions.push_back(params.aion(species.formula()));
        bions.push_back(params.bion(species.formula()));
    }

    // Collect the Debye-Huckel parameter b of the neutral species
    for(Index i : ineutral_species)
    {
        const auto species = mixture.species(i);
        bneutral.push_back(params.bneutral(species.formula()));
    }

    // The state of the aqueous mixture
    AqueousMixtureState state;

    // Define the activity model function of the aqueous mixture
    ActivityPropsFn fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x, extra] = args;

        // Evaluate the state of the aqueous mixture
        state = mixture.state(T, P, x);

        // Export the aqueous mixture and its state via the extra argument
        extra = { mixture, state };

        // Auxiliary constant references
        const auto& m = state.m;             // the molalities of all species
        const auto& ms = state.ms;           // the stoichiometric molalities of the charged species
        const auto& I = state.Is;            // the stoichiometric ionic strength
        const auto& rho = state.rho/1000;    // the density of water (in g/cm3)
        const auto& epsilon = state.epsilon; // the dielectric constant of water

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Auxiliary variables
		const auto ln_m = m.log();
		const auto xw = x[iwater];
		const auto ln_xw = log(xw);
		const auto mSigma = nwo * (1 - xw)/xw;
		const auto I2 = I*I;
		const auto sqrtI = sqrt(I);
		const auto sqrt_rho = sqrt(rho);
		const auto T_epsilon = T * epsilon;
		const auto sqrt_T_epsilon = sqrt(T_epsilon);
		const auto A = 1.824829238e+6 * sqrt_rho/(T_epsilon*sqrt_T_epsilon);
		const auto B = 50.29158649 * sqrt_rho/sqrt_T_epsilon;
		const auto sigmacoeff = (2.0/3.0)*A*I*sqrtI;

        // Set the first contribution to the activity of water
        ln_a[iwater] = mSigma;

        // Loop over all charged species in the aqueous mixture
        for(Index i = 0; i < num_charged_species; ++i)
        {
            // The index of the current charged species
            const auto ispecies = icharged_species[i];

            // The stoichiometric molality of the charged species
            const auto msi = ms[i];

            // The electrical charge of the charged species
            const auto z = charges[i];

            // Update the Lambda parameter of the Debye-Huckel activity coefficient model
            const auto Lambda = 1.0 + aions[i]*B*sqrtI;

			// Update the sigma parameter of the current ion
            const real sigma = (aions[i] != 0.0) ? 3.0*pow(Lambda - 1, -3) * ((Lambda - 1)*(Lambda - 3) + 2*log(Lambda)) : 2.0;

            // Calculate the ln activity coefficient of the current charged species
            ln_g[ispecies] = ln10 * (-A*z*z*sqrtI/Lambda + bions[i]*I);

            // Calculate the ln activity of the current charged species
            ln_a[ispecies] = ln_g[ispecies] + ln_m[ispecies];

            // Calculate the contribution of current ion to the ln activity of water
			ln_a[iwater] += msi*ln_g[ispecies] + sigmacoeff*sigma*ln10 - I2*bions[i]/(z*z)*ln10;
        }

        // Finalize the computation of the activity of water (in mole fraction scale)
        ln_a[iwater] *= -1.0/nwo;

        // Set the activity coefficient of water (mole fraction scale)
        ln_g[iwater] = ln_a[iwater] - ln_xw;

        // Loop over all neutral species in the aqueous mixture
        for(Index i = 0; i < num_neutral_species; ++i)
        {
            // The index of the current neutral species
            const auto ispecies = ineutral_species[i];

            // Calculate the ln activity coefficient of the current neutral species
            ln_g[ispecies] = ln10 * bneutral[i] * I;

            // Calculate the ln activity coefficient of the current neutral species
            ln_a[ispecies] = ln_g[ispecies] + ln_m[ispecies];
        }
    };

    return fn;
}

} // namespace detail

auto ActivityModelDebyeHuckelParams::aion(const ChemicalFormula& ion) const -> real
{
    return detail::get(ion, aions, aiondefault);
}

auto ActivityModelDebyeHuckelParams::bion(const ChemicalFormula& ion) const -> real
{
    return detail::get(ion, bions, biondefault);
}

auto ActivityModelDebyeHuckelParams::bneutral(const ChemicalFormula& neutral) const -> real
{
    return detail::get(neutral, bneutrals, bneutraldefault);
}

auto ActivityModelDebyeHuckelParams::setLimitingLaw() -> void
{
    aiondefault = 0.0;
    biondefault = 0.0;
    aions = {};
    bions = {};
}

auto ActivityModelDebyeHuckelParams::setKielland() -> void
{
    aiondefault = 0.0;
    biondefault = 0.0;
    aions = detail::aions_kielland;
    bions = {};
}

auto ActivityModelDebyeHuckelParams::setWATEQ4F() -> void
{
    aions = detail::aions_wateq4f;
    bions = detail::bions_wateq4f;
}

auto ActivityModelDebyeHuckelParams::setPHREEQC() -> void
{
    aions = detail::aions_phreeqc;
    bions = detail::bions_phreeqc;
    bneutraldefault = 0.1;
}

auto ActivityModelDebyeHuckel() -> ActivityModel
{
    return ActivityModelDebyeHuckelPHREEQC();
}

auto ActivityModelDebyeHuckel(ActivityModelDebyeHuckelParams params) -> ActivityModel
{
    return [=](const SpeciesList& species)
    {
        return detail::activityPropsFnDebyeHuckel(species, params);
    };
}

auto ActivityModelDebyeHuckelLimitingLaw() -> ActivityModel
{
    ActivityModelDebyeHuckelParams params;
    params.setLimitingLaw();
    return ActivityModelDebyeHuckel(params);
}

auto ActivityModelDebyeHuckelKielland() -> ActivityModel
{
    ActivityModelDebyeHuckelParams params;
    params.setKielland();
    return ActivityModelDebyeHuckel(params);
}

auto ActivityModelDebyeHuckelPHREEQC() -> ActivityModel
{
    ActivityModelDebyeHuckelParams params;
    params.setPHREEQC();
    return ActivityModelDebyeHuckel(params);
}

auto ActivityModelDebyeHuckelWATEQ4F() -> ActivityModel
{
    ActivityModelDebyeHuckelParams params;
    params.setWATEQ4F();
    return ActivityModelDebyeHuckel(params);
}

} // namespace Reaktoro
