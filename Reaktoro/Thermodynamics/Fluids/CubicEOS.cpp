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

#include "CubicEOS.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>
#include <Reaktoro/Math/Roots.hpp>

//=================================================================================================
// == REFERENCE ==
//=================================================================================================
// The implementation of this module is based on the presentation of the textbook:
//
//     Introduction to Chemical Engineering Thermodynamics, 8th Edition, 2017,
//     J.M. Smith, H. Van Ness, M. Abbott, M. Swihart
//
// More specifically, it is based on the information of the following chapters:
//
//      3.6 CUBIC EQUATIONS OF STATE, page 95
//     13.6 RESIDUAL PROPERTIES BY CUBIC EQUATIONS OF STATE, page 87
//
//-------------------------------------------------------------------------------------------------
// For more details and derivation of some formulas, check the document `notes/CubicEOS.lyx`.
//=================================================================================================

namespace Reaktoro {

using std::abs;
using std::log;
using std::sqrt;

const auto R = universalGasConstant;

namespace detail {

using AlphaResult = Tuple<real, real, real>;

/// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT for a given EOS.
auto alpha(CubicEOSModel type) -> Fn<AlphaResult(real, real, real)>
{
    // The alpha function for van der Waals EOS (see Table 3.1 of Smith et al. 2017)
    auto alphaVDW = [](real Tr, real TrT, real omega) -> AlphaResult
    {
        const real alpha = 1.0;
        const real alphaT = 0.0;
        const real alphaTT = 0.0;
        return { alpha, alphaT, alphaTT };
    };

    // The alpha function for Redlich-Kwong EOS
    auto alphaRK = [](real Tr, real TrT, real omega) -> AlphaResult
    {
        const real alpha = 1.0/sqrt(Tr);
        const real alphaTr = -0.5/Tr * alpha;
        const real alphaTrTr = -0.5/Tr * (alphaTr - alpha/Tr);
        const real alphaT = alphaTr*TrT;
        const real alphaTT = alphaTrTr*TrT*TrT;
        return { alpha, alphaT, alphaTT };
    };

    // The alpha function for Soave-Redlich-Kwong EOS
    auto alphaSRK = [](real Tr, real TrT, real omega) -> AlphaResult
    {
        const real m = 0.480 + 1.574*omega - 0.176*omega*omega;
        const real sqrtTr = sqrt(Tr);
        const real aux = 1.0 + m*(1.0 - sqrtTr);
        const real auxTr = -0.5*m/sqrtTr;
        const real auxTrTr = 0.25*m/(Tr*sqrtTr);
        const real alpha = aux*aux;
        const real alphaTr = 2.0*aux*auxTr;
        const real alphaTrTr = 2.0*(auxTr*auxTr + aux*auxTrTr);
        const real alphaT = alphaTr * TrT;
        const real alphaTT = alphaTrTr * TrT*TrT;
        return { alpha, alphaT, alphaTT };
    };

    // The alpha function for Peng-Robinson (1978) EOS
    auto alphaPR = [](real Tr, real TrT, real omega) -> AlphaResult
    {
        // Jaubert, J.-N., Vitu, S., Mutelet, F. and Corriou, J.-P., 2005.
        // Extension of the PPR78 model (predictive 1978, Peng–Robinson EOS
        // with temperature dependent kij calculated through a group
        // contribution method) to systems containing aromatic compounds.
        // Fluid Phase Equilibria, 237(1-2), pp.193–211.
        const real m = omega <= 0.491 ?
            0.374640 + 1.54226*omega - 0.269920*omega*omega :
            0.379642 + 1.48503*omega - 0.164423*omega*omega + 0.016666*omega*omega*omega;
        const real sqrtTr = sqrt(Tr);
        const real aux = 1.0 + m*(1.0 - sqrtTr);
        const real auxTr = -0.5*m/sqrtTr;
        const real auxTrTr = 0.25*m/(Tr*sqrtTr);
        const real alpha = aux*aux;
        const real alphaTr = 2.0*aux*auxTr;
        const real alphaTrTr = 2.0*(auxTr*auxTr + aux*auxTrTr);
        const real alphaT = alphaTr * TrT;
        const real alphaTT = alphaTrTr * TrT*TrT;
        return { alpha, alphaT, alphaTT };
    };

    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return alphaVDW;
        case CubicEOSModel::RedlichKwong: return alphaRK;
        case CubicEOSModel::SoaveRedlichKwong: return alphaSRK;
        case CubicEOSModel::PengRobinson: return alphaPR;
        default: return alphaPR;
    }
}

auto sigma(CubicEOSModel type) -> real
{
    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return 0.0;
        case CubicEOSModel::RedlichKwong: return 1.0;
        case CubicEOSModel::SoaveRedlichKwong: return 1.0;
        case CubicEOSModel::PengRobinson: return 1.0 + 1.4142135623730951;
        default: return 1.0 + 1.4142135623730951;
    }
}

auto epsilon(CubicEOSModel type) -> real
{
    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return 0.0;
        case CubicEOSModel::RedlichKwong: return 0.0;
        case CubicEOSModel::SoaveRedlichKwong: return 0.0;
        case CubicEOSModel::PengRobinson: return 1.0 - 1.4142135623730951;
        default: return 1.0 - 1.4142135623730951;
    }
}

auto Omega(CubicEOSModel type) -> real
{
    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return 1.0/8.0;
        case CubicEOSModel::RedlichKwong: return 0.08664;
        case CubicEOSModel::SoaveRedlichKwong: return 0.08664;
        case CubicEOSModel::PengRobinson: return 0.0777960739;
        default: return 0.0777960739;
    }
}

auto Psi(CubicEOSModel type) -> real
{
    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return 27.0/64.0;
        case CubicEOSModel::RedlichKwong: return 0.42748;
        case CubicEOSModel::SoaveRedlichKwong: return 0.42748;
        case CubicEOSModel::PengRobinson: return 0.457235529;
        default: return 0.457235529;
    }
}

/// Compute the local minimum of pressure along an isotherm of a cubic equation of state.
/// @param a The @eq{a_\mathrm{mix}} variable in the equation of state
/// @param b The @eq{b_\mathrm{mix}} variable in the equation of state
/// @param e The @eq{\epsilon} parameter in the cubic equation of state
/// @param s The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)
/// @return double
auto computeLocalMinimumPressureAlongIsotherm(double a, double b, double e, double s, double T) -> double
{
    const auto RT = R*T;

    auto V = b;
    auto Pprev = 0.0;

    const auto maxiters = 100;
    const auto tolerance = 1e-6;

    auto i = 0;
    for(; i < maxiters; ++i)
    {
        const auto t  = (V + e*b)*(V + s*b);
        const auto tV = 2*V + b*(e + s);

        const auto aux   = 1/b * sqrt(RT/(a*tV));
        const auto auxV  = -aux/tV;
        const auto auxVV = -3*auxV/tV;

        const auto q   = 1 + t*aux - V/b;
        const auto qV  = tV*aux + t*auxV - 1/b;
        const auto qVV = t*auxVV;

        const auto f = q*qV;
        const auto J = qV*qV + q*qVV;

        const auto dV = -f/J;

        V += dV;

        const auto P = RT/(V - b) - a/((V + e*b)*(V + s*b));

        if(abs(P - Pprev) < abs(P) * tolerance)
            return abs(q) < abs(qV) ? P : NaN;

        Pprev = P;
    }

    errorif(i == maxiters, "Could not compute the minimum pressure along an isotherm of a cubic equation of state.");

    return NaN;
}

/// Compute the residual Gibbs energy of the fluid for a given compressibility factor.
/// @param Z The compressibility factor
/// @param beta The @eq{\beta=Pb/(RT)} variable in the cubic equation of state
/// @param q The @eq{q=a/(bRT)} variable in the cubic equation of state
/// @param epsilon The @eq{\epsilon} parameter in the cubic equation of state
/// @param sigma The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)
auto computeResidualGibbsEnergy(double Z, double beta, double q, double epsilon, double sigma, double T) -> double
{
    auto I = 0.0;

    if(epsilon != sigma) // CASE I:  Eq. (13.72) of Smith et al. (2017)
        I = log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I=\frac{1}{\sigma-\epsilon}\ln\left(\frac{Z+\sigma\beta}{Z+\epsilon\beta}\right) }
    else // CASE II: Eq. (13.74) of Smith et al. (2017)
        I = beta/(Z + epsilon*beta); // @eq{ I=\frac{\beta}{Z+\epsilon\beta} }

    const auto Gres = R*T*(Z - 1 - log(Z - beta) - q*I); // from Eq. (13.74) of Smith et al. (2017)

    return Gres;
}

/// Determine the state of matter of the fluid when three real roots are available (either liquid or gas).
/// @param Zmin The compressibility factor with minimum value
/// @param Zmax The compressibility factor with maximum value
/// @param beta The @eq{\beta=Pb/(RT)} variable in the cubic equation of state
/// @param q The @eq{q=a/(bRT)} variable in the cubic equation of state
/// @param epsilon The @eq{\epsilon} parameter in the cubic equation of state
/// @param sigma The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)
/// @return StateOfMatter The state of matter of the fluid, by comparing the residual Gibbs energy of the two states.
auto determinePhysicalStateThreeRealRoots(double Zmin, double Zmax, double beta, double q, double epsilon, double sigma, double T) -> StateOfMatter
{
    const auto Gresmin = computeResidualGibbsEnergy(Zmin, beta, q, epsilon, sigma, T);
    const auto Gresmax = computeResidualGibbsEnergy(Zmax, beta, q, epsilon, sigma, T);
    return Gresmin < Gresmax ? StateOfMatter::Liquid : StateOfMatter::Gas;
}

/// Determine the state of matter of the fluid when only one real root is available (either supercritical or low pressure gas).
/// @param a The @eq{a_\mathrm{mix}} variable in the equation of state
/// @param b The @eq{b_\mathrm{mix}} variable in the equation of state
/// @param e The @eq{\epsilon} parameter in the cubic equation of state
/// @param s The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)
/// @param P The pressure (in Pa)
/// @return StateOfMatter The state of matter of the fluid, by comparing the residual Gibbs energy of the two states.
auto determinePhysicalStateOneRealRoot(double a, double b, double e, double s, double T, double P) -> StateOfMatter
{
    const auto Pmin = computeLocalMinimumPressureAlongIsotherm(a, b, e, s, T);
    return (Pmin != Pmin) ? StateOfMatter::Supercritical : (P < Pmin) ? StateOfMatter::Gas : StateOfMatter::Liquid;
}

} // namespace detail

struct CubicEOS::Impl
{
    /// The number of species in the phase.
    unsigned nspecies;

    /// The type of the cubic equation of state.
    CubicEOSModel model = CubicEOSModel::PengRobinson;

    /// The critical temperatures of the substances (in K).
    ArrayXr Tcr;

    /// The critical pressures of the substances (in Pa).
    ArrayXr Pcr;

    /// The acentric factor of the substances.
    ArrayXr omega;

    /// The function that calculates the interaction parameters kij and its temperature derivatives.
    CubicEOSInteractionParamsFn interaction_params_fn;

    ArrayXr a;     ///< Auxiliary array
    ArrayXr aT;    ///< Auxiliary array
    ArrayXr aTT;   ///< Auxiliary array
    ArrayXr b;     ///< Auxiliary array
    ArrayXr abar;  ///< Auxiliary array
    ArrayXr abarT; ///< Auxiliary array
    ArrayXr bbar;  ///< Auxiliary array

    /// Construct a CubicEOS::Impl instance.
    Impl(const Args& args)
    : nspecies(args.nspecies),
      Tcr(args.Tcr),
      Pcr(args.Pcr),
      omega(args.omega),
      model(args.model),
      interaction_params_fn(args.interaction_params_fn),
      a(args.nspecies),
      aT(args.nspecies),
      aTT(args.nspecies),
      b(args.nspecies),
      abar(args.nspecies),
      abarT(args.nspecies),
      bbar(args.nspecies)
    {

        assert(Tcr.size() == nspecies);
        assert(Pcr.size() == nspecies);
        assert(omega.size() == nspecies);
        assert(Tcr.minCoeff() > 0.0);
        assert(Pcr.minCoeff() > 0.0);
    }

    auto compute(CubicEOSProps& props, real T, real P, ArrayXrConstRef x) -> void
    {
        // Check if the mole fractions are zero or non-initialized
        if(x.size() == 0 || x.maxCoeff() <= 0.0)
            return;

        // Auxiliary variables
        const auto Psi = detail::Psi(model);
        const auto Omega = detail::Omega(model);
        const auto epsilon = detail::epsilon(model);
        const auto sigma = detail::sigma(model);
        const auto alphafn = detail::alpha(model);

        // Calculate the parameters `a` and `b` of the cubic equation of state for each species
        for(auto k = 0; k < nspecies; ++k)
        {
            const real factor = Psi*R*R*(Tcr[k]*Tcr[k])/Pcr[k]; // factor in Eq. (3.45) multiplying alpha
            const auto TrT = 1.0/Tcr[k];
            const auto Tr = T * TrT;
            const auto [alpha, alphaT, alphaTT] = alphafn(Tr, TrT, omega[k]);
            a[k]   = factor*alpha; // see Eq. (3.45)
            aT[k]  = factor*alphaT;
            aTT[k] = factor*alphaTT;
            b[k]   = Omega*R*Tcr[k]/Pcr[k]; // Eq. (3.44)
        }

        // Calculate the binary interaction parameters and its temperature derivatives
        CubicEOSInteractionParams ip;
        if(interaction_params_fn)
            ip = interaction_params_fn({T, a, aT, aTT, b});

        // Calculate the parameter `amix` of the phase and the partial molar parameters `abar` of each species
        real amix = {};
        real amixT = {};
        real amixTT = {};
        abar.fill(0.0);
        abarT.fill(0.0);
        for(auto i = 0; i < nspecies; ++i)
        {
            for(auto j = 0; j < nspecies; ++j)
            {
                const real r   = ip.k.size()   ? 1.0 - ip.k(i, j)   : real(1.0);
                const real rT  = ip.kT.size()  ?     - ip.kT(i, j)  : real(0.0);
                const real rTT = ip.kTT.size() ?     - ip.kTT(i, j) : real(0.0);

                const real s   = sqrt(a[i]*a[j]); // Eq. (13.93)
                const real sT  = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
                const real sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;

                const real aij   = r*s;
                const real aijT  = rT*s + r*sT;
                const real aijTT = rTT*s + 2.0*rT*sT + r*sTT;

                amix   += x[i] * x[j] * aij; // Eq. (13.92) of Smith et al. (2017)
                amixT  += x[i] * x[j] * aijT;
                amixTT += x[i] * x[j] * aijTT;

                abar[i]  += 2 * x[j] * aij;  // see Eq. (13.94)
                abarT[i] += 2 * x[j] * aijT;
            }
        }

        // Finalize the calculation of `abar` and `abarT`
        for(auto i = 0; i < nspecies; ++i)
        {
            abar[i] -= amix;
            abarT[i] -= amixT;
        }

        // Calculate the parameters bba[i] and bmix of the cubic equation of state
        //     bbar[i] = Omega*R*Tc[i]/Pc[i] as shown in Eq. (3.44)
        //     bmix = sum(x[i] * bbar[i])
        real bmix = {};
        for(auto i = 0; i < nspecies; ++i)
        {
            bbar[i] = Omega*R*Tcr[i]/Pcr[i]; // see Eq. (13.95) and unnumbered equation before Eq. (13.99)
            bmix += x[i] * bbar[i];  // Eq. (13.91) of Smith et al. (2017)
        }

        // Calculate the temperature and pressure derivatives of bmix
        const auto bmixT = 0.0; // no temperature dependence!
        const auto bmixP = 0.0; // no pressure dependence!

        // Calculate the auxiliary parameter beta and its partial derivatives betaT (at const P) and betaP (at const T)
        const real beta = P*bmix/(R*T); // Eq. (3.46)
        const real betaT = -beta/T; // note bmixT = 0
        const real betaP =  beta/P; // note bmixP = 0

        // Compute the auxiliary variable q and its partial derivatives qT, qTT (at const P) and qP (at const T)
        const real q = amix/(bmix*R*T); // Eq. (3.47)
        const real qT = q*(amixT/amix - 1.0/T); // === amixT/(bmix*R*T) - amix/(bmix*R*T*T)
        const real qTT = qT*qT/q + q*(amixTT/amix - amixT*amixT/(amix*amix) + 1.0/(T*T)); // === qT*(amixT/amix - 1.0/T) + q*(amixTT/amix - amixT*amixT/amix/amix + 1.0/T/T)
        const real qP = 0.0; // from Eq. (3.47), (dq/dP)_T := 0

        // Convert Eq. (3.48) into a cubic polynomial Z^3 + AZ^2 + BZ + C = 0, and compute the coefficients A, B, C of the cubic equation of state
        const real A = (epsilon + sigma - 1)*beta - 1;
        const real B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon + sigma - q)*beta;
        const real C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

        // Calculate AT := (dA/dT)_P, BT := (dB/dT)_P and CT := (dC/dT)_P (i.e., partial derivatives of A, B, C with respect to T at constant P)
        const real AT = (epsilon + sigma - 1)*betaT;
        const real BT = (epsilon*sigma - epsilon - sigma)*(2*beta*betaT) + qT*beta - (epsilon + sigma - q)*betaT;
        const real CT = -epsilon*sigma*(3*beta*beta*betaT) - qT*beta*beta - (epsilon*sigma + q)*(2*beta*betaT);

        // Calculate AP := (dA/dP)_T, BP := (dB/dP)_T and CP := (dC/dP)_T (i.e., partial derivatives of A, B, C with respect to P at constant T)
        const real AP = (epsilon + sigma - 1)*betaP;
        const real BP = (epsilon*sigma - epsilon - sigma)*(2*beta*betaP) + qP*beta - (epsilon + sigma - q)*betaP;
        const real CP = -epsilon*sigma*(3*beta*beta*betaP) - qP*beta*beta - (epsilon*sigma + q)*(2*beta*betaP);

        // Calculate cubic roots using cardano's method
        auto roots = realRoots(cardano(A, B, C));

        // Ensure there are either 1 or 3 real roots!
        assert(roots.size() == 1 || roots.size() == 3);

        // Determine the physical state of the fluid phase for given TPx conditions and its compressibility factor
        real Z = {};

        if(roots.size() == 3)
        {
            const auto Zmax = std::max({roots[0], roots[1], roots[2]});
            const auto Zmin = std::min({roots[0], roots[1], roots[2]});
            props.som = detail::determinePhysicalStateThreeRealRoots(Zmin, Zmax, beta, q, epsilon, sigma, T);
            Z = (props.som == StateOfMatter::Gas) ? Zmax : Zmin;
        }
        else
        {
            props.som = detail::determinePhysicalStateOneRealRoot(amix, bmix, epsilon, sigma, T, P);
            Z = roots[0];
        }

        // Calculate ZT := (dZ/dT)_P and ZP := (dZ/dP)_T
        const real ZT = -(AT*Z*Z + BT*Z + CT)/(3*Z*Z + 2*A*Z + B); // === (ZZZ + A*ZZ + B*Z + C)_T = 3*ZZ*ZT + AT*ZZ + 2*A*Z*ZT + BT*Z + B*ZT + CT = 0 => (3*ZZ + 2*A*Z + B)*ZT = -(AT*ZZ + BT*Z + CT)
        const real ZP = -(AP*Z*Z + BP*Z + CP)/(3*Z*Z + 2*A*Z + B); // === (ZZZ + A*ZZ + B*Z + C)_P = 3*ZZ*ZP + AP*ZZ + 2*A*Z*ZP + BP*Z + B*ZP + CP = 0 => (3*ZZ + 2*A*Z + B)*ZP = -(AP*ZZ + BP*Z + CP)

        //=========================================================================================
        // Calculate the integration factor I, IT := (dI/dT)_P and IP := (dI/dP)_T
        //=========================================================================================
        real I = {};
        real IT = {};
        real IP = {};

        if(epsilon != sigma) // CASE I:  Eq. (13.72) of Smith et al. (2017)
        {
            I = log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I=\frac{1}{\sigma-\epsilon}\ln\left(\frac{Z+\sigma\beta}{Z+\epsilon\beta}\right) }
            IT = ((ZT + sigma*betaT)/(Z + sigma*beta) - (ZT + epsilon*betaT)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I_{T}\equiv\left(\frac{\partial I}{\partial T}\right)_{P}=\frac{1}{\sigma-\epsilon}\left(\frac{Z_{T}+\sigma\beta_{T}}{Z+\sigma\beta}-\frac{Z_{T}+\epsilon\beta_{T}}{Z+\epsilon\beta}\right) }
            IP = ((ZP + sigma*betaP)/(Z + sigma*beta) - (ZP + epsilon*betaP)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I_{P}\equiv\left(\frac{\partial I}{\partial P}\right)_{T}=\frac{1}{\sigma-\epsilon}\left(\frac{Z_{P}+\sigma\beta_{P}}{Z+\sigma\beta}-\frac{Z_{P}+\epsilon\beta_{P}}{Z+\epsilon\beta}\right) }
        }
        else // CASE II: Eq. (13.74) of Smith et al. (2017)
        {
            I = beta/(Z + epsilon*beta); // @eq{ I=\frac{\beta}{Z+\epsilon\beta} }
            IT = I*(betaT/beta - (ZT + epsilon*betaT)/(Z + epsilon*beta)); // @eq{ I_{T}\equiv\left(\frac{\partial I}{\partial T}\right)_{P}=I\left(\frac{\beta_{T}}{\beta}-\frac{Z_{T}+\epsilon\beta_{T}}{Z+\epsilon\beta}\right) }
            IP = I*(betaP/beta - (ZP + epsilon*betaP)/(Z + epsilon*beta)); // @eq{ I_{P}\equiv\left(\frac{\partial I}{\partial P}\right)_{T}=I\left(\frac{\beta_{P}}{\beta}-\frac{Z_{P}+\epsilon\beta_{P}}{Z+\epsilon\beta}\right) }
        }

        //=========================================================================================
        // Calculate the ideal volume properties of the phase
        //=========================================================================================
        const real V0  =  R*T/P;
        const real V0T =  V0/T;
        const real V0P = -V0/P;

        //=========================================================================================
        // Calculate the corrected volumetric properties of the phase
        //=========================================================================================
        const auto& V  = props.V  = Z*V0;
        const auto& VT = props.VT = ZT*V0 + Z*V0T;
        const auto& VP = props.VP = ZP*V0 + Z*V0P;

        //=========================================================================================
        // Calculate the residual properties of the phase
        //=========================================================================================
        const auto& Gres  = props.Gres  = R*T*(Z - 1 - log(Z - beta) - q*I); // from Eq. (13.74) of Smith et al. (2017)
        const auto& Hres  = props.Hres  = R*T*(Z - 1 + T*qT*I); // equation after Eq. (13.74), but using T*qT instead of Tr*qTr, which is equivalent
        const auto& Cpres = props.Cpres = Hres/T + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT); // from Eq. (2.19), Cp(res) := (dH(res)/dT)P === R*(Z - 1 + T*qT*I) + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT) = H_res/T + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT)

        //=========================================================================================
        // Calculate the fugacity coefficients for each species
        //=========================================================================================
        for(auto k = 0; k < nspecies; ++k)
        {
            const real betak = P*bbar[k]/(R*T);
            const real qk    = (1 + abar[k]/amix - bbar[k]/bmix)*q;
            const real Ak    = (epsilon + sigma - 1.0)*betak - 1.0;
            const real Bk    = ((epsilon*sigma - epsilon - sigma)*(2*betak - beta) + qk - q)*beta - (epsilon + sigma - q)*betak;
            const real Ck    = (epsilon*sigma*(2*beta + 1) + 2*q - qk)*beta*beta - (2*(epsilon*sigma + q) + 3*epsilon*sigma*beta)*beta*betak;
            const real Zk    = -(Ak*Z*Z + (B + Bk)*Z + 2*C + Ck)/(3*Z*Z + 2*A*Z + B);

            const real Ik = (epsilon != sigma) ?
                I + ((Zk + sigma*betak)/(Z + sigma*beta) - (Zk + epsilon*betak)/(Z + epsilon*beta))/(sigma - epsilon) :
                I * (1 + betak/beta - (Zk + epsilon*betak)/(Z + epsilon*beta));

            props.ln_phi[k] = Zk - (Zk - betak)/(Z - beta) - log(Z - beta) + q*I - qk*I - q*Ik;
        }
    }
};

CubicEOS::CubicEOS(const Args& args)
: pimpl(new Impl(args))
{
}

CubicEOS::CubicEOS(const CubicEOS& other)
: pimpl(new Impl(*other.pimpl))
{}

CubicEOS::~CubicEOS()
{}

auto CubicEOS::operator=(CubicEOS other) -> CubicEOS&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto CubicEOS::setModel(CubicEOSModel model) -> void
{
    pimpl->model = model;
}

auto CubicEOS::setInteractionParamsFunction(const CubicEOSInteractionParamsFn& func) -> void
{
    pimpl->interaction_params_fn = func;
}

auto CubicEOS::compute(CubicEOSProps& props, real T, real P, ArrayXrConstRef x) -> void
{
    return pimpl->compute(props, T, P, x);
}

} // namespace Reaktoro
