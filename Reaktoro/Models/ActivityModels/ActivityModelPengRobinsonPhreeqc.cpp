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
#include <Reaktoro/Singletons/CriticalProps.hpp>

namespace Reaktoro {
namespace {

using std::log;
using std::sqrt;

/// Returns true if string `str` starts with `sub0`, or `sub1`, etc.
template<typename... Subs>
auto startswith(String const& str, const char* sub0, Subs... subs)
{
    if constexpr(sizeof...(Subs) > 0)
        return str.find(sub0) == 0 || startswith(str, subs...);
    return str.find(sub0) == 0;
}

/// The arguments for the evaluation of the properties a fluid phase using the Peng-Robinson-Phreeqc equation of state model.
struct PengRobinsonPhreeqcArgs
{
    /// The temperature of the fluid phase (in K).
    real const& TK;

    /// The pressure of the fluid phase (in atm).
    real const& Patm;

    /// The mole fractions of the species in the fluid phase (in mol/mol).
    ArrayXrConstRef const& xgases;

    /// The critical temperatures of the species in the fluid phase (in K).
    ArrayXrConstRef const& Tcr;

    /// The critical pressures of the species in the fluid phase (in atm).
    ArrayXrConstRef const& Pcr;

    /// The acentric factors of the species in the fluid phase.
    ArrayXrConstRef const& omega;

    /// The names of the species in the fluid phase.
    Strings const& gases;
};

/// The properties computed for a fluid phase using the Peng-Robinson-Phreeqc equation of state model.
struct PengRobinsonPhreeqcProps
{
    /// The molar volume of the fluid phase (in liter/mol)
    real V_m;

    /// The fugacity coefficients of the fluid species (in natural log scale).
    ArrayXr ln_phi;
};

// Based on LDBLE calc_PR(std::vector<class phase*> phase_ptrs, LDBLE P, LDBLE TK, LDBLE V_m);
/* ---------------------------------------------------------------------- */
/*  Calculate fugacity and fugacity coefficient for gas pressures if critical T and P
    are defined.
  1) Solve molar volume V_m or total pressure P from Peng-Robinson's EOS:
  P = R * T / (V_m - b) - a * aa / (V_m^2 + 2 * b * V_m - b^2)
     a = 0.457235 * (R * T_c)^2 / P_c
     b = 0.077796 * R * T_c / P_c
     aa = (1 + kk * (1 - T_r^0.5))^2
     kk = 0.37464 + 1.54226 * omega - 0.26992 * omega^2
     T_r = T / T_c
  multicomponent gas phase:
     use: b_sum = Sum(x_i * b), x_i is mole-fraction
          a_aa_sum = Sum_i( Sum_j(x_i * x_j * (a_i * aa_i * a_j * aa_j)^0.5) )
  2) Find the fugacity coefficient phi for gas i:
  log(phi_i) = B_ratio * (z - 1) - log(z - B) + A / (2.8284 * B) * (B_ratio - 2 / a_aa_sum * a_aa_sum2) *\
           log((z + 2.4142 * B) / (z - 0.4142 * B))
     B_ratio = b_i / b_sum
     A = a_aa_sum * P / R_TK^2
     B = b_sum * P / R_TK
     a_aa_sum2 = Sum_j(x_j * (a_aa_i * a_aa_j)^0.5
  3) correct the solubility of gas i with:
  pr_si_f = log10(phi_i) -  Delta_V_i * (P - 1) / (2.303 * R * TK);
*/


void evaluatePengRobinsonPhreeqc(PengRobinsonPhreeqcArgs const& args, PengRobinsonPhreeqcProps& props)
{
    const auto R_LITER_ATM = 0.0820597; // L-atm/deg-mol
    int i, i1, n_g = args.gases.size();
    real T_c, P_c;
    real A, B, B_r, kk, oo, a_aa, T_r;
    real a_aa_sum2;
    real phi;
    real /*R_TK,*/ R = R_LITER_ATM; /* L atm / (K mol) */
    real r3[4], r3_12, rp, rp3, rq, rz, ri, ri1, one_3 = 0.33333333333333333;
    real disct, vinit, v1, ddp, dp_dv, dp_dv2;
    int it;
    // class phase *phase_ptr, *phase_ptr1;
    // cxxGasPhase * gas_phase_ptr = use.Get_gas_phase_ptr();
    bool halved;
    real P = args.Patm;
    real TK = args.TK;
    real R_TK = R * TK;
    real b_sum;
    real a_aa_sum;

    auto& V_m = props.V_m;
    auto& ln_phi = props.ln_phi;

    ArrayXr pr_a(n_g);
    ArrayXr pr_b(n_g);
    ArrayXr pr_alpha(n_g);
    ArrayXr pr_aa_sum2(n_g);

    for(i = 0; i < n_g; i++)
    {
        errorif(args.Tcr[i] == 0.0, "Expecting a non-zero critical temperature for gas ", args.gases[i], " in Peng-Robinson model.");
        errorif(args.Pcr[i] == 0.0, "Expecting a non-zero critical pressure for gas ", args.gases[i], " in Peng-Robinson model.");

        T_c = args.Tcr[i];
        P_c = args.Pcr[i];
        pr_a[i] = 0.457235 * R * R * T_c * T_c / P_c;
        pr_b[i] = 0.077796 * R * T_c / P_c;
        T_r = TK / T_c;
        oo = args.omega[i];
        kk = 0.37464 + oo * (1.54226 - 0.26992 * oo);
        pr_alpha[i] = pow(1 + kk * (1 - sqrt(T_r)), 2);
    }

    for (i = 0; i < n_g; i++)
    {
        a_aa_sum2 = 0.0;
        b_sum += args.xgases[i] * pr_b[i];

        for (i1 = 0; i1 < n_g; i1++)
        {
            if(args.xgases[i1] == 0.0)
                continue;

            a_aa = sqrt(pr_a[i]*pr_alpha[i] * pr_a[i1]*pr_alpha[i1]);

            if(startswith(args.gases[i], "H2O"))
            {
                if(startswith(args.gases[i1], "CO2"))
                    a_aa *= 0.81; // Soreide and Whitson, 1992, FPE 77, 217
                else if(startswith(args.gases[i1], "H2S"))
                    a_aa *= 0.81;
                else if(startswith(args.gases[i1], "CH4", "Mtg", "Methane", "METHANE"))
                    a_aa *= 0.51;
                else if(startswith(args.gases[i1], "N2", "Ntg"))
                    a_aa *= 0.51;
                else if(startswith(args.gases[i1], "C2H6", "Ethane", "ETHANE"))
                    a_aa *= 0.51;
                else if(startswith(args.gases[i1], "C3H8", "Propane", "PROPANE"))
                    a_aa *= 0.45;
            }
            if(startswith(args.gases[i1], "H2O"))
            {
                if(startswith(args.gases[i], "CO2"))
                    a_aa *= 0.81;
                else if(startswith(args.gases[i], "H2S"))
                    a_aa *= 0.81;
                else if(startswith(args.gases[i], "CH4", "Mtg", "Methane", "METHANE"))
                    a_aa *= 0.51;
                else if(startswith(args.gases[i], "N2", "Ntg"))
                    a_aa *= 0.51;
                else if(startswith(args.gases[i], "C2H6", "Ethane", "ETHANE"))
                    a_aa *= 0.51;
                else if(startswith(args.gases[i], "C3H8", "Propane", "PROPANE"))
                    a_aa *= 0.45;
            }
            a_aa_sum += args.xgases[i] * args.xgases[i1] * a_aa;
            a_aa_sum2 += args.xgases[i1] * a_aa;
        }
        pr_aa_sum2[i] = a_aa_sum2;
    }

    real b2 = b_sum * b_sum;

    if(P < 1e-10)
        P = 1e-10;

    r3[1] = b_sum - R_TK / P;
    r3_12 = r3[1] * r3[1];
    r3[2] = -3.0 * b2 + (a_aa_sum - R_TK * 2.0 * b_sum) / P;
    r3[3] = b2 * b_sum + (R_TK * b2 - b_sum * a_aa_sum) / P;
    // solve t^3 + rp*t + rq = 0.
    // molar volume V_m = t - r3[1] / 3...
    rp = r3[2] - r3_12 / 3;
    rp3 = rp * rp * rp;
    rq = (2.0 * r3_12 * r3[1] - 9.0 * r3[1] * r3[2]) / 27 + r3[3];
    rz = rq * rq / 4 + rp3 / 27;
    if(rz >= 0) // Cardono's method...
    {
        ri = sqrt(rz);
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
        ri = sqrt(- rp3 / 27); // rp < 0
        ri1 = acos(- rq / 2 / ri);
        V_m = 2.0 * pow(ri, one_3) * cos(ri1 / 3) - r3[1] / 3;
    }

    // calculate the fugacity coefficients...
    for (i = 0; i < n_g; i++)
    {
        if(args.xgases[i] == 0.0)
        {
            ln_phi[i] = 0.0;
            continue;
        }
        rz = P * V_m / R_TK;
        A = a_aa_sum * P / (R_TK * R_TK);
        B = b_sum * P / R_TK;
        B_r = pr_b[i] / b_sum;
        if(rz > B)
            ln_phi[i] = B_r * (rz - 1) - log(rz - B) + A / (2.828427 * B) * (B_r - 2.0 * pr_aa_sum2[i] / a_aa_sum) * log((rz + 2.41421356 * B) / (rz - 0.41421356 * B));
        else ln_phi[i] = -3.0; // fugacity coefficient = 0.05
    }
}

} // namespace anonymous

auto createActivityModelPengRobinsonPhreeqc(SpeciesList const& specieslist) -> ActivityModel
{
    // Ensure all species are gas, liquid, or fluid
    // for(auto species : specieslist)
    //     errorifnot(oneof(species.aggregateState(), AggregateState::Gas, AggregateState::Liquid, AggregateState::Fluid), "Expecting all species for Peng-Robinson-Phreeqc model to have a gas, liquid, or fluid aggregate state but this is not the case for ", species.name(), ".");

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

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& [T, P, x] = args;

        const auto Patm = P / atmToPascal; // convert from Pa to atm
        const auto Pbar = P / barToPascal; // convert from Pa to bar

        PengRobinsonPhreeqcArgs prargs{ T, Patm, x, Tcr, Pcr, omega, names };

        evaluatePengRobinsonPhreeqc(prargs, prprops);

        props.Vx   = prprops.V_m * 1e-3; // convert from liter/mol to m3/mol
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
