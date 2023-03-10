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

#include "PhreeqcWater.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Memoization.hpp>

namespace Reaktoro {
namespace PhreeqcUtils {

auto waterThermoProps(real T, real P) -> PhreeqcWaterThermoProps
{
    //--------------------------------------------------------------------------------
    // Implementation based on PHREEQC method `Phreeqc::calc_rho_0`
    //--------------------------------------------------------------------------------

    // Compute density of pure water using Wagner and Pruss, 2002, JPCRD 31,
    // 387, eqn. 2.6, along the saturation pressure line + interpolation 0 -
    // 300 oC, 0.006 - 1000 atm... just like in PHREEQC.

    PhreeqcWaterThermoProps props;

    auto& [rho_0, kappa_0] = props; // references to the properties computed below

    const auto pascal_to_atm = 9.86923e-6;

    auto tc = T - 273.15;
    auto pa = P * pascal_to_atm;

    warning(tc > 350.0, "The water density model used in PHREEQC is "
        "valid up to 350 degC, but given temperature is ", tc, " degC! "
        "Setting temperature to 350 degC, just like PHREEQC does.");

    if(tc > 350.0)
        tc = 350.0;

    // Eq. (2.6) of Wagner and Pruss (2002)
    const auto Tc = 647.096;
    const auto th = 1 - T / Tc;
    const auto b1 = 1.99274064;
    const auto b2 = 1.09965342;
    const auto b3 = -0.510839303;
    const auto b4 = -1.75493479;
    const auto b5 = -45.5170352;
    const auto b6 = -6.7469445e5;

    // Saturated liquid density of water (in kg/m3)
    const auto rho_0_sat = 322.0 * (1.0 + b1 * pow(th, 1./3.) + b2 * pow(th, 2./3.) + b3 * pow(th, 5./3.) +
               b4 * pow(th, 16./3.) + b5 * pow(th, 43./3.) + b6 * pow(th, 110./3));

    // Pressure interpolation
    const auto p0 =  5.1880000E-02 + tc * (-4.1885519E-04 + tc * ( 6.6780748E-06 + tc * (-3.6648699E-08 + tc *  8.3501912E-11)));
    const auto p1 = -6.0251348E-06 + tc * ( 3.6696407E-07 + tc * (-9.2056269E-09 + tc * ( 6.7024182E-11 + tc * -1.5947241E-13)));
    const auto p2 = -2.2983596E-09 + tc * (-4.0133819E-10 + tc * ( 1.2619821E-11 + tc * (-9.8952363E-14 + tc *  2.3363281E-16)));
    const auto p3 =  7.0517647E-11 + tc * ( 6.8566831E-12 + tc * (-2.2829750E-13 + tc * ( 1.8113313E-15 + tc * -4.2475324E-18)));

    // Note: PHREEQC uses activity of water as a factor in the saturation pressure.
    // Not sure why. This is not used here.
    const auto p_sat = exp(11.6702 - 3816.44 / (T - 46.13));

    // Note: PHREEQC sets pressure to saturation pressure if given pressure is
    // below saturation. This is most likely done because at this condition,
    // water is in vapor state, but here a liquid density is needed.
    if(pa < p_sat)
        pa = p_sat;

    pa -= (p_sat - 1e-6);

    // The density of water at T,P (in kg/m3)
    rho_0 = rho_0_sat + pa*(p0 + pa*(p1 + pa*(p2 + sqrt(pa)*p3)));

    if(rho_0 < 0.01)
        rho_0 = 0.01;

    // Compute the water compressibility, d(ln(rho))/dP (in 1/atm)
    kappa_0 = (p0 + pa * (2*p1 + pa*(3*p2 + sqrt(pa)*3.5*p3))) / rho_0;

    // Convert water density from kg/m3 to g/cm3 (like in PHREEQC)
    rho_0 /= 1e3;

    return props;
}

auto waterElectroProps(real T, real P, PhreeqcWaterThermoProps wtp) -> PhreeqcWaterElectroProps
{
    // COMMENT FROM PHREEQC:
    // Relative dielectric constant of pure water, eps as a function of (P, T)
    // Bradley and Pitzer, 1979, JPC 83, 1599.
    // (newer data in Fernandez et al., 1995, JPCRD 24, 33,
    // and Fernandez et al., 1997, JPCRD 26, 1125, show its correctness)
    // + d(eps)/d(P), Debye-Hueckel A and B, and Av (for Av, see Pitzer et al., 1984, JPCRD 13, p. 4)

    PhreeqcWaterElectroProps props;

    auto& [eps_r, DH_A, DH_B, DH_Av, ZBrn, QBrn] = props; // references to the electro properties computed below

    const auto& [rho_0, kappa_0] = wtp; // references to the given thermo properties of water

    const auto pascal_to_atm = 9.86923e-6;

    auto tc = T - 273.15;
    auto pa = P * pascal_to_atm;

    if(tc > 350.0)
        tc = 350.0;

    const auto u1 =  3.4279e2;
    const auto u2 = -5.0866e-3;
    const auto u3 =  9.469e-7;
    const auto u4 = -2.0525;
    const auto u5 =  3.1159e3;
    const auto u6 = -1.8289e2;
    const auto u7 = -8.0325e3;
    const auto u8 =  4.2142e6;
    const auto u9 =  2.1417;

    const auto d1000 = u1 * exp(T * (u2 + T * u3)); // relative dielectric constant at 1000 bar

    const auto c = u4 + u5 / (u6 + T);
    const auto b = u7 + u8 / T + u9 * T;
    const auto pb = pa * 1.01325; // pa in bar

    eps_r = d1000 + c * log((b + pb) / (b + 1e3)); // relative dielectric constant

    if(eps_r <= 0.0)
    {
        eps_r = 10.0;
        warning(true, "Relative dielectric constant is negative. "
            "Temperature is out of range of parameterization. "
            "Setting it to 10.0 just like PHREEQC.");
    }

    /* qe^2 / (eps_r * kB * T) = 4.803204e-10**2 / 1.38065e-16 / (eps_r * T)
                               = 1.671008e-3 (esu^2 / (erg/K)) / (eps_r * T) */
    const auto e2_DkT = 1.671008e-3 / (eps_r * T);
    const auto pi = 3.14159265358979; // taken from PHREEQC
    const auto LOG_10 = ln10;
    const auto AVOGADRO = 6.02252e23; // taken from PHREEQC
    const auto R_LITER_ATM = 0.0820597; // taken from PHREEQC

    DH_B = sqrt(8 * pi * AVOGADRO * e2_DkT * rho_0 / 1e3);  // Debye length parameter, 1/cm(mol/kg)^-0.5

    DH_A = DH_B * e2_DkT / (2. * LOG_10); //(mol/kg)^-0.5

    /* Debye-Hueckel limiting slope = DH_B *  e2_DkT * RT * (d(ln(eps_r)) / d(P) - compressibility) */
    DH_Av = DH_B * e2_DkT * R_LITER_ATM * 1e3 * T * (c / (b + pb) * 1.01325 / eps_r - kappa_0 / 3.); // (cm3/mol)(mol/kg)^-0.5

    DH_B /= 1e8; // kappa, 1/Angstrom(mol/kg)^-0.5

    /* the Born functions, * 41.84 to give molal volumes in cm3/mol... */
    ZBrn = (- 1 / eps_r + 1.0) * 41.84004;
    QBrn = c / (b + pb) / eps_r / eps_r * 41.84004;

    return props;
}

auto waterProps(real T, real P) -> PhreeqcWaterProps
{
    const auto wtp = waterThermoProps(T, P);
    const auto wep = waterElectroProps(T, P, wtp);
    return {wtp, wep};
}

auto waterPropsMemoized(real T, real P) -> PhreeqcWaterProps
{
    static thread_local auto memoized_water_props = memoizeLast(waterProps);
    return memoized_water_props(T, P);
}

auto waterDensityPhreeqc(real T, real P) -> real
{
    const auto Tmax = 350.0 + 273.15;

    warning(T > Tmax, "Evaluating water density for Debye-Huckel coefficients at temperature ", T, " K, but maximum supported is ", Tmax, " K.");

    // Critical properties of water as given in eqs. (2.2a-2.2c) in Wagner and Pruss (2002)
    const auto Tcr = 647.096;  // critical temperature (in K)
    const auto Pcr = 22.064e6; // critical presssure (in Pa)
    const auto Dcr = 322.0;    // critical density (in kg/m3)

    // The theta variable in Wagner and Pruss (2002) defined below eq. (2.5a)
    const auto theta = 1 - T / Tcr;

    // The coefficients b1 to b6 in eq. (2.6) of Wagner and Pruss (2002)
    const auto b1 =  1.99274064;
    const auto b2 =  1.09965342;
    const auto b3 = -0.510839303;
    const auto b4 = -1.75493479;
    const auto b5 = -45.5170352;
    const auto b6 = -6.7469445e5;

    // Calculate the density of water along the liquid saturation curve
    const auto Dsat = Dcr * (1.0 +
            b1 * pow(theta, 1./3.) +
            b2 * pow(theta, 2./3.) +
            b3 * pow(theta, 5./3.) +
            b4 * pow(theta, 16./3.) +
            b5 * pow(theta, 43./3.) +
            b6 * pow(theta, 110./3));

    // Factor to convert from Pa to atm
    const auto pascal_to_atm = 9.86923e-6;

    // Calculate temperature in C and and presssure in atm
    const auto Tc = T - 273.15;
    const auto Pa = P * pascal_to_atm;

    // Interpolate the density of water from saturation curve to desired pressure
    // This has been taken from PHREEQC, LDBLE Phreeqc::calc_rho_0(LDBLE tc, LDBLE pa) method in utilities.cpp
    const auto p0 =  5.1880000E-02 + Tc * (-4.1885519E-04 + Tc * ( 6.6780748E-06 + Tc * (-3.6648699E-08 + Tc *  8.3501912E-11)));
    const auto p1 = -6.0251348E-06 + Tc * ( 3.6696407E-07 + Tc * (-9.2056269E-09 + Tc * ( 6.7024182E-11 + Tc * -1.5947241E-13)));
    const auto p2 = -2.2983596E-09 + Tc * (-4.0133819E-10 + Tc * ( 1.2619821E-11 + Tc * (-9.8952363E-14 + Tc *  2.3363281E-16)));
    const auto p3 =  7.0517647E-11 + Tc * ( 6.8566831E-12 + Tc * (-2.2829750E-13 + Tc * ( 1.8113313E-15 + Tc * -4.2475324E-18)));

    // Calculate the pressure along the liquid saturation curve corresponding to given temperature
    const auto Psat = exp(11.6702 - 3816.44 / (T - 46.13)); // in atm

    // Determine the pressure value in effect, which is either given P (in atm) or saturated pressure (the minimum pressure to maintain liquid state)
    const auto Peff = max(Pa, Psat);

    // Calculate delta(P) = Peff - P(sat) (all in atm)
    const auto dP = Peff - Psat + 1e-6; // add small pressure value to avoid zero

    // Calculate the interpolation of water density from saturation value to a corrected one for given pressure
    const auto D = Dsat + dP * (p0 + dP * (p1 + dP * (p2 + sqrt(dP) * p3)));

    return D;
}

auto waterDielectricConstantPhreeqc(real T, real P) -> real
{
    const auto Tmax = 350.0 + 273.15;

    warning(T > Tmax, "Evaluating water dielectric constant for Debye-Huckel coefficients at temperature ", T, " K, but maximum supported is ", Tmax, " K.");

    const auto u1 =  3.4279e2;
    const auto u2 = -5.0866e-3;
    const auto u3 =  9.469e-7;
    const auto u4 = -2.0525;
    const auto u5 =  3.1159e3;
    const auto u6 = -1.8289e2;
    const auto u7 = -8.0325e3;
    const auto u8 =  4.2142e6;
    const auto u9 =  2.1417;

    const auto epsilon1000 = u1 * exp(T * (u2 + T * u3)); // relative dielectric constant at 1000 bar
    const auto c = u4 + u5 / (u6 + T);
    const auto b = u7 + u8 / T + u9 * T;
    const auto pb = P * 1e-5; // pressure from Pa to bar

    const auto epsilon = epsilon1000 + c * log((b + pb) / (b + 1e3)); // relative dielectric constant

    return epsilon;
}

} // namespace PhreeqcUtils
} // namespace Reaktoro
