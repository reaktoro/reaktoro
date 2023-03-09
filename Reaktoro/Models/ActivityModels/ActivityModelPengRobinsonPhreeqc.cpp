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

#include "ActivityModelPengRobinsonPhreeqc.hpp"

// C++ includes
#include <cassert>
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Singletons/CriticalProps.hpp>

namespace Reaktoro {
namespace {

using std::log;
using std::sqrt;

/// The arguments for the evaluation of the properties a fluid phase using the Peng-Robinson-Phreeqc equation of state model.
struct PengRobinsonPhreeqcArgs
{
    /// The temperature of the fluid phase (in K).
    real const& TK;

    /// The pressure of the fluid phase (in atm).
    real const& Patm;

    /// The mole fractions of the species in the fluid phase (in mol/mol).
    ArrayXrConstRef const& x;
};

/// The properties computed for a fluid phase using the Peng-Robinson-Phreeqc equation of state model.
struct PengRobinsonPhreeqcProps
{
    /// The molar volume of the fluid phase (in liter/mol)
    real V_m;

    /// The fugacity coefficients of the fluid species (in natural log scale).
    ArrayXr ln_phi;
};

/// The arguments for the evaluation of the properties a fluid phase using the Peng-Robinson-Phreeqc equation of state model.
struct PengRobinsonPhreeqc
{
    /// The names of the species in the fluid phase.
    const Strings species;

    /// The critical temperatures of the species in the fluid phase (in K).
    const ArrayXr Tcr;

    /// The critical pressures of the species in the fluid phase (in atm).
    const ArrayXr Pcr;

    /// The acentric factors of the species in the fluid phase.
    const ArrayXr omega;

    // Auxiliary working arrays

    ArrayXr pr_a;
    ArrayXr pr_b;
    ArrayXr pr_alpha;
    ArrayXr pr_aa_sum2;

    /// Evaluates the Peng-Robinson model similar to that implemented in PHREEQC.
    auto evaluate(PengRobinsonPhreeqcArgs const& args, PengRobinsonPhreeqcProps& props)
    {
        const auto R_LITER_ATM = 0.0820597; // L-atm/deg-mol
        const auto R = R_LITER_ATM; // L atm / (K mol)
        const auto one_3 = 1.0/3.0;
        const auto n_g = species.size();

        real P = args.Patm;
        real TK = args.TK;
        real R_TK = R * TK;

        real b_sum;
        real a_aa_sum;
        real a_aa_sum2;

        pr_a.resize(n_g);
        pr_b.resize(n_g);
        pr_alpha.resize(n_g);
        pr_aa_sum2.resize(n_g);

        auto& V_m = props.V_m;
        auto& ln_phi = props.ln_phi;

        for(auto i = 0; i < n_g; ++i)
        {
            const auto T_c = Tcr[i];
            const auto P_c = Pcr[i];
            const auto oo = omega[i];
            const auto T_r = TK / T_c;
            const auto kk = 0.37464 + oo * (1.54226 - 0.26992 * oo);
            pr_a[i] = 0.457235 * R * R * T_c * T_c / P_c;
            pr_b[i] = 0.077796 * R * T_c / P_c;
            pr_alpha[i] = pow(1 + kk * (1 - sqrt(T_r)), 2);
        }

        for(auto i = 0; i < n_g; ++i)
        {
            a_aa_sum2 = 0.0;
            b_sum += args.x[i] * pr_b[i];

            for(auto i1 = 0; i1 < n_g; ++i1)
            {
                if(args.x[i1] == 0.0)
                    continue;

                auto a_aa = sqrt(pr_a[i]*pr_alpha[i] * pr_a[i1]*pr_alpha[i1]);

                if(startswith(species[i], "H2O"))
                {
                    if(startswith(species[i1], "CO2"))
                        a_aa *= 0.81; // Soreide and Whitson, 1992, FPE 77, 217
                    else if(startswith(species[i1], "H2S"))
                        a_aa *= 0.81;
                    else if(startswith(species[i1], "CH4", "Mtg", "Methane", "METHANE"))
                        a_aa *= 0.51;
                    else if(startswith(species[i1], "N2", "Ntg"))
                        a_aa *= 0.51;
                    else if(startswith(species[i1], "C2H6", "Ethane", "ETHANE"))
                        a_aa *= 0.51;
                    else if(startswith(species[i1], "C3H8", "Propane", "PROPANE"))
                        a_aa *= 0.45;
                }
                if(startswith(species[i1], "H2O"))
                {
                    if(startswith(species[i], "CO2"))
                        a_aa *= 0.81;
                    else if(startswith(species[i], "H2S"))
                        a_aa *= 0.81;
                    else if(startswith(species[i], "CH4", "Mtg", "Methane", "METHANE"))
                        a_aa *= 0.51;
                    else if(startswith(species[i], "N2", "Ntg"))
                        a_aa *= 0.51;
                    else if(startswith(species[i], "C2H6", "Ethane", "ETHANE"))
                        a_aa *= 0.51;
                    else if(startswith(species[i], "C3H8", "Propane", "PROPANE"))
                        a_aa *= 0.45;
                }
                a_aa_sum += args.x[i] * args.x[i1] * a_aa;
                a_aa_sum2 += args.x[i1] * a_aa;
            }
            pr_aa_sum2[i] = a_aa_sum2;
        }

        real b2 = b_sum * b_sum;

        if(P < 1e-10)
            P = 1e-10;

        real r3[4];

        r3[1] = b_sum - R_TK / P;
        r3[2] = -3.0 * b2 + (a_aa_sum - R_TK * 2.0 * b_sum) / P;
        r3[3] = b2 * b_sum + (R_TK * b2 - b_sum * a_aa_sum) / P;

        // solve t^3 + rp*t + rq = 0.
        // molar volume V_m = t - r3[1] / 3...
        const auto r3_12 = r3[1] * r3[1];
        const auto rp = r3[2] - r3_12 / 3;
        const auto rp3 = rp * rp * rp;
        const auto rq = (2.0 * r3_12 * r3[1] - 9.0 * r3[1] * r3[2]) / 27 + r3[3];
        const auto rz = rq * rq / 4 + rp3 / 27;

        if(rz >= 0) // Cardono's method...
        {
            auto ri = sqrt(rz);
            if(ri + rq / 2 <= 0)
            {
                V_m = pow(ri - rq / 2, one_3) + pow(- ri - rq / 2, one_3) - r3[1] / 3;
            }
            else
            {
                ri = - pow(ri + rq / 2, one_3);
                V_m = ri - rp / (3.0 * ri) - r3[1] / 3;
            }
        }
        else // use complex plane...
        {
            const auto ri = sqrt(- rp3 / 27); // rp < 0
            const auto ri1 = acos(- rq / 2 / ri);
            V_m = 2.0 * pow(ri, one_3) * cos(ri1 / 3) - r3[1] / 3;
        }

        // calculate the fugacity coefficients...
        for(auto i = 0; i < n_g; ++i)
        {
            if(args.x[i] == 0.0)
            {
                ln_phi[i] = 0.0;
                continue;
            }
            const auto rz = P * V_m / R_TK;
            const auto A = a_aa_sum * P / (R_TK * R_TK);
            const auto B = b_sum * P / R_TK;
            const auto B_r = pr_b[i] / b_sum;
            if(rz > B)
                ln_phi[i] = B_r * (rz - 1) - log(rz - B) + A / (2.828427 * B) * (B_r - 2.0 * pr_aa_sum2[i] / a_aa_sum) * log((rz + 2.41421356 * B) / (rz - 0.41421356 * B));
            else ln_phi[i] = -3.0; // fugacity coefficient = 0.05
        }
    }
};

} // namespace anonymous

auto createActivityModelPengRobinsonPhreeqc(SpeciesList const& specieslist) -> ActivityModel
{
    const auto numspecies = specieslist.size();

    Strings names(numspecies); // the names of the species
    ArrayXr Tcr(numspecies);   // the critical temperatures of the species (in K)
    ArrayXr Pcr(numspecies);   // the critical pressures of the species (in atm)
    ArrayXr omega(numspecies); // the acentric factors of the species

    for(auto [i, species] : enumerate(specieslist))
    {
        const auto crprops = CriticalProps::get({
            species.substance(),
            species.formula(),
            species.name()
        });

        errorifnot(crprops.has_value(),
            "Could not find critical properties for substance ", species.formula().str(), " "
            "while creating a Peng-Robinson-Phreeqc model. ",
            "The following are ways to fix this error:\n"
            "  - use CriticalProps::append(crprops) to register the critical properties of this substance in the CriticalProps database,\n"
            "  - use CriticalProps::setMissingAs(\"He\") to consider all missing substances in the CriticalProps database as if they were He.");

        names[i] = species.name();
        Tcr[i] = crprops->temperature();
        Pcr[i] = crprops->pressure() / atmToPascal; // convert from Pa to atm
        omega[i] = crprops->acentricFactor();
    }

    PengRobinsonPhreeqcProps prprops;
    prprops.ln_phi.resize(numspecies);

    PengRobinsonPhreeqc prphreeqc{ names, Tcr, Pcr, omega };

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& [T, P, x] = args;

        const auto Patm = P / atmToPascal; // convert from Pa to atm
        const auto Pbar = P / barToPascal; // convert from Pa to bar

        prphreeqc.evaluate({ T, Patm, x}, prprops);

        props.Vx   = prprops.V_m * 0.001; // convert from liter/mol to m3/mol
        props.VxT  = 0.0; // PHREEQC does not compute this property
        props.VxP  = 0.0; // PHREEQC does not compute this property
        props.Gx   = 0.0; // PHREEQC does not compute this property
        props.Hx   = 0.0; // PHREEQC does not compute this property
        props.Cpx  = 0.0; // PHREEQC does not compute this property
        props.ln_g = prprops.ln_phi;
        props.ln_a = prprops.ln_phi + log(x) + log(Pbar);
        props.som  = StateOfMatter::Gas;
    };

    return fn;
}

auto ActivityModelPengRobinsonPhreeqc() -> ActivityModelGenerator
{
    return [=](SpeciesList const& species) { return createActivityModelPengRobinsonPhreeqc(species); };
}

} // namespace Reaktoro
