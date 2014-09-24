// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#include "AqueousElectroStateUtils.hpp"

// C++ includes
#include <cmath>

#include "WaterThermoStateUtils.hpp"
using std::pow;
using std::log;
using std::abs;

// Reaktor includes
#include <Reaktor/Common/ConvertUtils.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Thermodynamics/AqueousElectroState.hpp>
#include <Reaktor/Thermodynamics/WaterThermoState.hpp>
#include <Reaktor/Thermodynamics/WaterThermoStateUtils.hpp>

namespace Reaktor {
namespace internal {

/// The \eta constant in the HKF model (in units of (A*cal)/mol)
const double eta = 1.66027e+05;

} /* namespace internal */

using namespace internal;

auto functionG(double T, double P, const WaterThermoState& wts) -> FunctionG
{
    // The function G
    FunctionG funcG;

    // The temperature in units of celsius and pressure in units of bar
    const double TdegC = convert<degK,degC>(T);
    const double Pbar  = convert<Pa,bar>(P);

    // Check if the point (T,P) is inside region III or the shaded region in Fig. 6 of
    // Shock and others (1992), on page 809. In this case, we assume the g function to be zero.
    if(wts.density > 1000.0 or wts.density < 350.0)
        return funcG;

    // Auxiliary references
    auto& g   = funcG.g;
    auto& gT  = funcG.gT;
    auto& gP  = funcG.gP;
    auto& gTT = funcG.gTT;
    auto& gTP = funcG.gTP;
    auto& gPP = funcG.gPP;

    // Use equations (24)-(31) of Shock and others (1992) to compute `g` and its derivatives on region I
    const double ag1 = -2.037662;
    const double ag2 =  5.747000e-03;
    const double ag3 = -6.557892e-06;

    const double bg1 =  6.107361;
    const double bg2 = -1.074377e-02;
    const double bg3 =  1.268348e-05;

    const double ag = ag1 + ag2*TdegC + ag3*TdegC*TdegC;
    const double bg = bg1 + bg2*TdegC + bg3*TdegC*TdegC;

    const double agT = ag2 + 2*ag3*TdegC;
    const double bgT = bg2 + 2*bg3*TdegC;

    const double agTT = 2*ag3;
    const double bgTT = 2*bg3;

    const double r =  wts.density/1000.0;

    const double alpha  = -wts.densityT/wts.density;
    const double beta   =  wts.densityP/wts.density;
    const double alphaT = -wts.densityTT/wts.density + alpha*alpha;
    const double alphaP = -wts.densityTP/wts.density - alpha*beta;
    const double betaP  =  wts.densityPP/wts.density - beta*beta;

    g   =  ag * pow(1 - r, bg);
    gT  =   g * (agT/ag + bgT*log(1 - r) + r*alpha*bg/(1 - r));
    gP  =  -g * r*beta*bg/(1 - r);
    gTT =   g * (agTT/ag - pow(agT/ag, 2) + bgTT*log(1 - r) + r*alpha*bg/(1 - r) * (2*bgT/bg + alphaT/alpha - alpha - r*alpha/(1 - r))) + gT*gT/g;
    gTP =  gP * (bgT/bg - alpha - alphaP/beta - r*alpha/(1 - r)) + gP*gT/g;
    gPP =  gP * (gP/g + beta + betaP/beta + r*beta/(1 - r));

    // Check if the point (T,P) is inside region II, as depicted in Fig. 6 of Shock and others (1992), on page 809
    if(TdegC > 155.0 and TdegC < 355.0 and Pbar < 1000.0)
    {
        // Use equations (32)-(44) of Shock and others (1992) to compute the function g and its partial derivatives on region II
        const double af1 =  3.666660e+01; // unit: K
        const double af2 = -1.504956e-10; // unit: A:bar^{-3}
        const double af3 =  5.017990e-14; // unit: A:bar^{-4}

        const double auxT  = (TdegC - 155)/300;
        const double auxT1 = pow(auxT, 4.8);
        const double auxT2 = pow(auxT, 16);

        const double auxP  = 1000 - Pbar;
        const double auxP1 = pow(auxP, 3);
        const double auxP2 = pow(auxP, 4);

        const double ft   = auxT1 + af1*auxT2;
        const double ftT  = ( 4.80 * auxT1 +  16.0 * af1*auxT2)/(300 * auxT);
        const double ftTT = (18.24 * auxT1 + 240.0 * af1*auxT2)/pow(300 * auxT, 2);

        const double fp   =  af2*auxP1 + af3*auxP2;
        const double fpP  = -(3.0 * af2*auxP1 +  4.0 * af3*auxP2)/(auxP*1.e+5);       // convert derivative from bar to Pa
        const double fpPP =  (6.0 * af2*auxP1 + 12.0 * af3*auxP2)/pow(auxP*1.e+5, 2); // convert derivative from bar^2 to Pa^2

        g   -= ft  * fp;
        gT  -= fp  * ftT;
        gP  -= ft  * fpP;
        gTT -= fp  * ftTT;
        gTP -= ftT * fpP;
        gPP -= ft  * fpPP;
    }

    return funcG;
}

auto aqueousEletroStateHKF(const FunctionG& g, const AqueousSpecies& species) -> AqueousElectroState
{
    // Get the HKF thermodynamic parameters of the aqueous species
    const auto& hkf = species.thermoparams.hkf();

    // The species electro instance to be calculated
    AqueousElectroState se;

    // Check if the aqueous species is neutral or the ion H+ and set the electrostatic data accordingly
    if(species.charge == 0.0 or species.name == "H+")
    {
        se.w   = hkf.wref;
        se.wT  = 0.0;
        se.wP  = 0.0;
        se.wTT = 0.0;
        se.wTP = 0.0;
        se.wPP = 0.0;
    }
    else
    {
        const double z    = species.charge;
        const double wref = hkf.wref;

        const double reref = z*z/(wref/eta + z/3.082);
        const double re    = reref + abs(z) * g.g;

        const double X1 =  -eta * (abs(z*z*z)/(re*re) - z/pow(3.082 + g.g, 2));
        const double X2 = 2*eta * (z*z*z*z/(re*re*re) - z/pow(3.082 + g.g, 3));

        se.re    = re;
        se.reref = reref;
        se.w     = eta * (z*z/re - z/(3.082 + g.g));
        se.wT    = X1 * g.gT;
        se.wP    = X1 * g.gP;
        se.wTT   = X1 * g.gTT + X2 * g.gT * g.gT;
        se.wTP   = X1 * g.gTP + X2 * g.gT * g.gP;
        se.wPP   = X1 * g.gPP + X2 * g.gP * g.gP;
    }

    return se;
}

auto aqueousEletroStateHKF(double T, double P, const AqueousSpecies& species) -> AqueousElectroState
{
    WaterThermoState wt = waterThermoStateWagnerPruss(T, P);

    FunctionG g = functionG(T, P, wt);

    return aqueousEletroStateHKF(g, species);
}

} // namespace Reaktor
