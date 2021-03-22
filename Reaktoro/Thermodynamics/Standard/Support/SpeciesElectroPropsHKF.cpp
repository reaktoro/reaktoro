// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "SpeciesElectroPropsHKF.hpp"

// C++ includes
#include <cmath>
using std::abs;
using std::log;
using std::pow;

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Thermodynamics/Standard/StandardThermoModelHKF.hpp>
#include <Reaktoro/Thermodynamics/Standard/Support/SpeciesElectroProps.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>

namespace Reaktoro {
namespace {

/// The @eq{\eta} constant in the HKF model (in units of A*(J/mol))
const double eta = 6.94656968e+05; // from HKF, 1.66027e+05 A*(cal/mol)

} // namespace

auto gHKF::compute(real T, real P, const WaterThermoState& wts) -> gHKF
{
    gHKF res;

    auto& [g, gT, gP, gTT, gTP, gPP] = res;

    // The temperature in units of celsius and pressure in units of bar
    const auto TdegC = T - 273.15;
    const auto Pbar  = P * 1.0e-5;

    // Check if the point (T,P) is inside region III or the shaded region in Fig. 6 of
    // Shock and others (1992), on page 809. In this case, we assume the g function to be zero.
    if(wts.density > 1000.0 || wts.density < 350.0)
        return res;

    // Use equations (24)-(31) of Shock and others (1992) to compute `g` and its derivatives on region I
    const auto ag1 = -2.037662;
    const auto ag2 =  5.747000e-03;
    const auto ag3 = -6.557892e-06;

    const auto bg1 =  6.107361;
    const auto bg2 = -1.074377e-02;
    const auto bg3 =  1.268348e-05;

    const auto ag = ag1 + ag2*TdegC + ag3*TdegC*TdegC;
    const auto bg = bg1 + bg2*TdegC + bg3*TdegC*TdegC;

    const auto agT = ag2 + 2*ag3*TdegC;
    const auto bgT = bg2 + 2*bg3*TdegC;

    const auto agTT = 2*ag3;
    const auto bgTT = 2*bg3;

    const auto r =  wts.density/1000.0;

    const auto alpha  = -wts.densityT/wts.density;
    const auto beta   =  wts.densityP/wts.density;
    const auto alphaT = -wts.densityTT/wts.density + alpha*alpha;
    const auto alphaP = -wts.densityTP/wts.density - alpha*beta;
    const auto betaP  =  wts.densityPP/wts.density - beta*beta;

    g   =  ag * pow(1 - r, bg);
    gT  =   g * (agT/ag + bgT*log(1 - r) + r*alpha*bg/(1 - r));
    gP  =  -g * r*beta*bg/(1 - r);
    gTT =   g * (agTT/ag - pow(agT/ag, 2) + bgTT*log(1 - r) + r*alpha*bg/(1 - r) * (2*bgT/bg + alphaT/alpha - alpha - r*alpha/(1 - r))) + gT*gT/g;
    gTP =  gP * (bgT/bg - alpha - alphaP/beta - r*alpha/(1 - r)) + gP*gT/g;
    gPP =  gP * (gP/g + beta + betaP/beta + r*beta/(1 - r));

    // Check if the point (T,P) is inside region II, as depicted in Fig. 6 of Shock and others (1992), on page 809
    if(TdegC > 155.0 && TdegC < 355.0 && Pbar < 1000.0)
    {
        // Use equations (32)-(44) of Shock and others (1992) to compute the function g and its partial derivatives on region II
        const auto af1 =  3.666660e+01; // unit: K
        const auto af2 = -1.504956e-10; // unit: A:bar^{-3}
        const auto af3 =  5.017990e-14; // unit: A:bar^{-4}

        const auto auxT  = (TdegC - 155)/300;
        const auto auxT1 = pow(auxT, 4.8);
        const auto auxT2 = pow(auxT, 16);

        const auto auxP  = 1000 - Pbar;
        const auto auxP1 = pow(auxP, 3);
        const auto auxP2 = pow(auxP, 4);

        const auto ft   = auxT1 + af1*auxT2;
        const auto ftT  = ( 4.80 * auxT1 +  16.0 * af1*auxT2)/(300 * auxT);
        const auto ftTT = (18.24 * auxT1 + 240.0 * af1*auxT2)/pow(300 * auxT, 2);

        const auto fp   =  af2*auxP1 + af3*auxP2;
        const auto fpP  = -(3.0 * af2*auxP1 +  4.0 * af3*auxP2)/(auxP*1.e+5);       // convert derivative from bar to Pa
        const auto fpPP =  (6.0 * af2*auxP1 + 12.0 * af3*auxP2)/pow(auxP*1.e+5, 2); // convert derivative from bar^2 to Pa^2

        g   -= ft  * fp;
        gT  -= fp  * ftT;
        gP  -= ft  * fpP;
        gTT -= fp  * ftTT;
        gTP -= ftT * fpP;
        gPP -= ft  * fpPP;
    }

    return res;
}

auto speciesElectroPropsHKF(const gHKF& gstate, const StandardThermoModelParamsHKF& params) -> SpeciesElectroProps
{
    // The species electro instance to be calculated
    SpeciesElectroProps se;

    // Auxiliary references
    const auto& [g, gT, gP, gTT, gTP, gPP] = gstate;

    // Check if the aqueous species is neutral or H+, and set its electrostatic data accordingly
    if(params.charge == 0.0 || params.formula.equivalent("H+"))
    {
        se.w   = params.wref;
        se.wT  = 0.0;
        se.wP  = 0.0;
        se.wTT = 0.0;
        se.wTP = 0.0;
        se.wPP = 0.0;
    }
    else
    {
        const auto z = params.charge;
        const auto wref = params.wref;

        const auto reref = z*z/(wref/eta + z/3.082);
        const auto re    = reref + abs(z) * g;

        const auto X1 =  -eta * (abs(z*z*z)/(re*re) - z/pow(3.082 + g, 2));
        const auto X2 = 2*eta * (z*z*z*z/(re*re*re) - z/pow(3.082 + g, 3));

        se.re    = re;
        se.reref = reref;
        se.w     = eta * (z*z/re - z/(3.082 + g));
        se.wT    = X1 * gT;
        se.wP    = X1 * gP;
        se.wTT   = X1 * gTT + X2 * gT * gT;
        se.wTP   = X1 * gTP + X2 * gT * gP;
        se.wPP   = X1 * gPP + X2 * gP * gP;
    }

    return se;
}

} // namespace Reaktoro
