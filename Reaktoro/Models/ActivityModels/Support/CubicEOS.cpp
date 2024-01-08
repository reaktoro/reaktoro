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
namespace CubicEOS {

using std::abs;
using std::log;
using std::sqrt;

const auto R = universalGasConstant;

// For implementations below, see Table 3.1 of Smith et al. (2017).

auto EquationModelVanDerWaals() -> EquationModel
{
    EquationModel eqmodel;
    eqmodel.sigma   = 0.0;
    eqmodel.epsilon = 0.0;
    eqmodel.Omega   = 1.0/8.0;
    eqmodel.Psi     = 27.0/64.0;

    eqmodel.alphafn = [](AlphaModelArgs const& args) -> Alpha
    {
        return { 1.0, 0.0, 0.0 };
    };

    return eqmodel;
}

auto EquationModelRedlichKwong() -> EquationModel
{
    EquationModel eqmodel;
    eqmodel.sigma   = 1.0;
    eqmodel.epsilon = 0.0;
    eqmodel.Omega   = 0.08664;
    eqmodel.Psi     = 0.42748;

    eqmodel.alphafn = [](AlphaModelArgs const& args) -> Alpha
    {
        auto const& [Tr, TrT, omega] = args;
        auto const alpha = 1.0/sqrt(Tr);
        auto const alphaTr = -0.5/Tr * alpha;
        auto const alphaTrTr = -0.5/Tr * (alphaTr - alpha/Tr);
        auto const alphaT = alphaTr*TrT;
        auto const alphaTT = alphaTrTr*TrT*TrT;
        return { alpha, alphaT, alphaTT };
    };

    return eqmodel;
}

auto EquationModelSoaveRedlichKwong() -> EquationModel
{
    EquationModel eqmodel;
    eqmodel.sigma   = 1.0;
    eqmodel.epsilon = 0.0;
    eqmodel.Omega   = 0.08664;
    eqmodel.Psi     = 0.42748;

    eqmodel.alphafn = [](AlphaModelArgs const& args) -> Alpha
    {
        auto const& [Tr, TrT, omega] = args;
        auto const m = 0.480 + 1.574*omega - 0.176*omega*omega;
        auto const sqrtTr = sqrt(Tr);
        auto const aux = 1.0 + m*(1.0 - sqrtTr);
        auto const auxTr = -0.5*m/sqrtTr;
        auto const auxTrTr = 0.25*m/(Tr*sqrtTr);
        auto const alpha = aux*aux;
        auto const alphaTr = 2.0*aux*auxTr;
        auto const alphaTrTr = 2.0*(auxTr*auxTr + aux*auxTrTr);
        auto const alphaT = alphaTr * TrT;
        auto const alphaTT = alphaTrTr * TrT*TrT;
        return { alpha, alphaT, alphaTT };
    };

    return eqmodel;
}

auto EquationModelPengRobinsonAux(Index year) -> EquationModel
{
    EquationModel eqmodel;
    eqmodel.sigma   = 1.0 + 1.4142135623730951;
    eqmodel.epsilon = 1.0 - 1.4142135623730951;
    eqmodel.Omega   = 0.0777960739;
    eqmodel.Psi     = 0.457235529;

    Fn<real(real const&)> mPR76 = [](real const& omega) -> real
    {
        return 0.374640 + 1.54226*omega - 0.269920*omega*omega;
    };

    Fn<real(real const&)> mPR78 = [](real const& omega) -> real
    {
        return omega < 0.491 ?
            0.374640 + 1.54226*omega - 0.269920*omega*omega :
            0.379642 + 1.48503*omega - 0.164423*omega*omega + 0.016666*omega*omega*omega;  // Jaubert, J.-N., Vitu, S., Mutelet, F. and Corriou, J.-P., 2005. Extension of the PPR78 model (predictive 1978, Peng–Robinson EOS with temperature dependent kij calculated through a group contribution method) to systems containing aromatic compounds. Fluid Phase Equilibria, 237(1-2), pp.193–211.
    };

    auto mPR = year == 76 ? mPR76 : mPR78; // Note: MSVC 19.29.30147.0 fails to compile with mPR76 and mPR78 are not explicitly assigned with type Fn<real(real const&)>!

    eqmodel.alphafn = [=](AlphaModelArgs const& args) -> Alpha
    {
        auto const& [Tr, TrT, omega] = args;
        auto const m = mPR(omega);
        auto const sqrtTr = sqrt(Tr);
        auto const aux = 1.0 + m*(1.0 - sqrtTr);
        auto const auxTr = -0.5*m/sqrtTr;
        auto const auxTrTr = 0.25*m/(Tr*sqrtTr);
        auto const alpha = aux*aux;
        auto const alphaTr = 2.0*aux*auxTr;
        auto const alphaTrTr = 2.0*(auxTr*auxTr + aux*auxTrTr);
        auto const alphaT = alphaTr * TrT;
        auto const alphaTT = alphaTrTr * TrT*TrT;
        return { alpha, alphaT, alphaTT };
    };

    return eqmodel;
}

auto EquationModelPengRobinson76() -> EquationModel
{
    return EquationModelPengRobinsonAux(76);
}

auto EquationModelPengRobinson78() -> EquationModel
{
    return EquationModelPengRobinsonAux(78);
}

auto EquationModelPengRobinson() -> EquationModel
{
    return EquationModelPengRobinson78();
}

namespace detail {

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

        const auto w   = 1/b * sqrt(RT/(a*tV));
        const auto wV  = -w/tV;
        const auto wVV = -3*wV/tV;

        const auto q   = 1 + t*w - V/b;
        const auto qV  = tV*w + t*wV - 1/b;
        const auto qVV = t*wVV;

        const auto f = q*qV;
        const auto J = qV*qV + q*qVV;

        const auto dV = -f/J;

        if(V == b && dV < 0.0)
            return NaN; // if V = b and dV is negative, consider supercritical (no local minimum for V > b)

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

/// Return the critical temperatures of the substances as an array, checking if their values are valid.
auto getCriticalTemperatures(Vec<Substance> const& substances) -> ArrayXr
{
    ArrayXr res(substances.size());
    auto i = 0; for(auto const& subs : substances) {
        errorif(std::isnan(subs.Tcr), "The critical temperature of `", subs.formula, "` was not explicitly initialized in CubicEOS::EquationSpecs::substances.");
        errorif(subs.Tcr <= 0.0, "The critical temperature of `", subs.formula, "` should be a positive number but it was set to ", subs.Tcr, " K.");
        res[i++] = subs.Tcr;
    }
    return res;
}

/// Return the critical pressures of the substances as an array, checking if their values are valid.
auto getCriticalPressures(Vec<Substance> const& substances) -> ArrayXr
{
    ArrayXr res(substances.size());
    auto i = 0; for(auto const& subs : substances) {
        errorif(std::isnan(subs.Tcr), "The critical pressure of `", subs.formula, "` was not explicitly initialized in CubicEOS::EquationSpecs::substances.");
        errorif(subs.Tcr <= 0.0, "The critical pressure of `", subs.formula, "` should be a positive number but it was set to ", subs.Pcr, " Pa.");
        res[i++] = subs.Pcr;
    }
    return res;
}

/// Return the acentric factors of the substances as an array, checking if their values are valid.
auto getAccentricFactors(Vec<Substance> const& substances) -> ArrayXr
{
    ArrayXr res(substances.size());
    auto i = 0; for(auto const& subs : substances) {
        errorif(std::isnan(subs.Tcr), "The acentric factor of `", subs.formula, "` was not explicitly initialized in CubicEOS::EquationSpecs::substances.");
        res[i++] = subs.omega;
    }
    return res;
}

/// Return the chemical formulas of the substances in the fluid phase.
auto createSubstanceList(Vec<Substance> const& substances) -> Strings
{
    return vectorize(substances, RKT_LAMBDA(x, x.formula));
}

} // namespace detail

struct Equation::Impl
{
    /// The specifications for the cubic equation of state.
    EquationSpecs const eqspecs;

    /// The number of species in the fluid phase.
    Index const nspecies;

    /// The critical temperatures of the species in the fluid phase (in K).
    ArrayXr const Tcr;

    /// The critical pressures of the species in the fluid phase (in Pa).
    ArrayXr const Pcr;

    /// The acentric factors of the species in the fluid phase.
    ArrayXr const omega;

    /// The chemical formulas of the substances in the fluid phase.
    Strings const substances;

    // Auxiliary arrays

    ArrayXr a;
    ArrayXr aT;
    ArrayXr aTT;
    ArrayXr alpha;
    ArrayXr alphaT;
    ArrayXr alphaTT;
    ArrayXr b;
    ArrayXr abar;
    ArrayXr abarT;
    ArrayXr bbar;
    Bip bip;

    /// Construct an Equation::Impl object.
    Impl(EquationSpecs const& eqspecs)
    : eqspecs(eqspecs),
      nspecies(eqspecs.substances.size()),
      Tcr(detail::getCriticalTemperatures(eqspecs.substances)),
      Pcr(detail::getCriticalPressures(eqspecs.substances)),
      omega(detail::getAccentricFactors(eqspecs.substances)),
      substances(detail::createSubstanceList(eqspecs.substances))
    {
        errorifnot(eqspecs.eqmodel.alphafn.initialized(), "The alpha function in CubicEOS::EquationSpecs::alphafn has not been initialized.");

        a       = zeros(nspecies);
        aT      = zeros(nspecies);
        aTT     = zeros(nspecies);
        alpha   = zeros(nspecies);
        alphaT  = zeros(nspecies);
        alphaTT = zeros(nspecies);
        b       = zeros(nspecies);
        abar    = zeros(nspecies);
        abarT   = zeros(nspecies);
        bbar    = zeros(nspecies);
        bip.k   = zeros(nspecies, nspecies);
        bip.kT  = zeros(nspecies, nspecies);
        bip.kTT = zeros(nspecies, nspecies);
    }

    auto compute(Props& props, real const& T, real const& P, ArrayXrConstRef const& x) -> void
    {
        // Check if the mole fractions are zero or non-initialized
        if(x.size() == 0 || x.maxCoeff() <= 0.0)
            return;

        // Auxiliary references
        auto const& sigma   = eqspecs.eqmodel.sigma.value();
        auto const& epsilon = eqspecs.eqmodel.epsilon.value();
        auto const& Omega   = eqspecs.eqmodel.Omega.value();
        auto const& Psi     = eqspecs.eqmodel.Psi.value();
        auto const& alphafn = eqspecs.eqmodel.alphafn;

        // Calculate the parameters `a` and `b` of the cubic equation of state for each species
        for(auto k = 0; k < nspecies; ++k)
        {
            const auto factor = Psi*R*R*(Tcr[k]*Tcr[k])/Pcr[k]; // factor in Eq. (3.45) multiplying alpha
            const auto TrT = 1.0/Tcr[k];
            const auto Tr = T * TrT;
            const auto [alphak, alphaTk, alphaTTk] = alphafn({ Tr, TrT, omega[k] });
            alpha[k]   = alphak;
            alphaT[k]  = alphaTk;
            alphaTT[k] = alphaTTk;
            a[k]       = factor*alphak; // see Eq. (3.45)
            aT[k]      = factor*alphaTk;
            aTT[k]     = factor*alphaTTk;
            b[k]       = Omega*R*Tcr[k]/Pcr[k]; // Eq. (3.44)
        }

        // Calculate the binary interaction parameters and its temperature derivatives
        if(eqspecs.bipmodel.initialized())
            eqspecs.bipmodel(bip, { substances, T, Tcr, Pcr, omega, a, aT, aTT, alpha, alphaT, alphaTT, b });

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
                auto const r   = 1.0 - bip.k(i, j);
                auto const rT  = -bip.kT(i, j);
                auto const rTT = -bip.kTT(i, j);

                auto const s   = sqrt(a[i]*a[j]); // Eq. (13.93)
                auto const sT  = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
                auto const sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;

                auto const aij   = r*s;
                auto const aijT  = rT*s + r*sT;
                auto const aijTT = rTT*s + 2.0*rT*sT + r*sTT;

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
        props.Vi.resize(nspecies);
        props.ln_phi.resize(nspecies);
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

            props.Vi[k] = R * T * Zk / P;
            props.ln_phi[k] = Zk - (Zk - betak)/(Z - beta) - log(Z - beta) + q*I - qk*I - q*Ik;
        }
    }
};

Equation::Equation(EquationSpecs const& eqspecs)
: pimpl(new Impl(eqspecs))
{
}

Equation::Equation(Equation const& other)
: pimpl(new Impl(*other.pimpl))
{}

Equation::~Equation()
{}

auto Equation::operator=(Equation other) -> Equation&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Equation::equationSpecs() const -> EquationSpecs const&
{
    return pimpl->eqspecs;
}

auto Equation::compute(Props& props, real const& T, real const& P, ArrayXrConstRef const& x) -> void
{
    return pimpl->compute(props, T, P, x);
}

auto BipModelPhreeqc(Strings const& substances, BipModelParamsPhreeqc const& params) -> BipModel
{
    auto isubstance = [&](auto... substrs)
    {
        return indexfn(substances, RKT_LAMBDA(substance, startswith(substance, substrs...)));
    };

    const auto iH2O  = isubstance("H2O");
    const auto iCO2  = isubstance("CO2");
    const auto iH2S  = isubstance("H2S");
    const auto iCH4  = isubstance("CH4", "Mtg", "Methane", "METHANE");
    const auto iN2   = isubstance("N2", "Ntg");
    const auto iC2H6 = isubstance("C2H6", "Ethane", "ETHANE");
    const auto iC3H8 = isubstance("C3H8", "Propane", "PROPANE");

    const auto size = substances.size();

    auto evalfn = [=](Bip& bip, BipModelArgs const& args) -> void
    {
        auto& [k, kT, kTT] = bip;

        if(iH2O < size)
        {
            if( iCO2 < size) k(iH2O,  iCO2) = k( iCO2, iH2O) = params.kH2O_CO2.value();
            if( iH2S < size) k(iH2O,  iH2S) = k( iH2S, iH2O) = params.kH2O_H2S.value();
            if( iCH4 < size) k(iH2O,  iCH4) = k( iCH4, iH2O) = params.kH2O_CH4.value();
            if(  iN2 < size) k(iH2O,   iN2) = k(  iN2, iH2O) = params.kH2O_N2.value();
            if(iC2H6 < size) k(iH2O, iC2H6) = k(iC2H6, iH2O) = params.kH2O_C2H6.value();
            if(iC3H8 < size) k(iH2O, iC3H8) = k(iC3H8, iH2O) = params.kH2O_C3H8.value();
        }
    };

    Data paramsdata;
    auto node = paramsdata["BipModel"]["Phreeqc"];
    node["kH2O_CO2"]  = params.kH2O_CO2;
    node["kH2O_H2S"]  = params.kH2O_H2S;
    node["kH2O_CH4"]  = params.kH2O_CH4;
    node["kH2O_N2"]   = params.kH2O_N2;
    node["kH2O_C2H6"] = params.kH2O_C2H6;
    node["kH2O_C3H8"] = params.kH2O_C3H8;

    return BipModel(evalfn, paramsdata);
}

auto BipModelSoreideWhitson(Strings const& substances, BipModelParamsSoreideWhitson const& params) -> BipModel
{
    auto isubstance = [&](auto... substrs)
    {
        return indexfn(substances, RKT_LAMBDA(substance, startswith(substance, substrs...)));
    };

    const auto iH2O    = isubstance("H2O");
    const auto iCO2    = isubstance("CO2");
    const auto iH2S    = isubstance("H2S");
    const auto iCH4    = isubstance("CH4", "Mtg", "Methane", "METHANE");
    const auto iN2     = isubstance("N2", "Ntg");
    const auto iC2H6   = isubstance("C2H6", "Ethane", "ETHANE");
    const auto iC3H8   = isubstance("C3H8", "Propane", "PROPANE");
    const auto inC4H10 = isubstance("nC4H10", "n-C4H10", "n-Butane", "n-Butane", "n-BUTANE", "N-BUTANE");

    const auto size = substances.size();

    auto evalfn = [=](Bip& bip, BipModelArgs const& args) -> void
    {
        auto& [k, kT, kTT] = bip;

        const auto TrH2S = args.T / args.Tcr[iH2S];

        if(iH2O < size)
        {
            if( iCO2   < size) k(iH2O,    iCO2) = k(   iCO2, iH2O) = params.kH2O_CO2.value();
            if(  iN2   < size) k(iH2O,     iN2) = k(    iN2, iH2O) = params.kH2O_N2.value();
            if( iH2S   < size) k(iH2O,    iH2S) = k(   iH2S, iH2O) = params.kH2O_H2S_a1.value() + params.kH2O_H2S_a2.value() * TrH2S;
            if( iCH4   < size) k(iH2O,    iCH4) = k(   iCH4, iH2O) = params.kH2O_CH4.value();
            if(iC2H6   < size) k(iH2O,   iC2H6) = k(  iC2H6, iH2O) = params.kH2O_C2H6.value();
            if(iC3H8   < size) k(iH2O,   iC3H8) = k(  iC3H8, iH2O) = params.kH2O_C3H8.value();
            if(inC4H10 < size) k(iH2O, inC4H10) = k(inC4H10, iH2O) = params.kH2O_nC4H10.value();
        }
    };

    Data paramsdata;
    auto node = paramsdata["BipModel"]["SoreideWhitson"];
    paramsdata["kH2O_CO2"] = params.kH2O_CO2;
    paramsdata["kH2O_N2"] = params.kH2O_N2;
    paramsdata["kH2O_CH4"] = params.kH2O_CH4;
    paramsdata["kH2O_C2H6"] = params.kH2O_C2H6;
    paramsdata["kH2O_C3H8"] = params.kH2O_C3H8;
    paramsdata["kH2O_nC4H10"] = params.kH2O_nC4H10;
    paramsdata["kH2O_H2S_a1"] = params.kH2O_H2S_a1;
    paramsdata["kH2O_H2S_a2"] = params.kH2O_H2S_a2;

    return BipModel(evalfn, paramsdata);
}

// DEPRECATED METHODS TO BE REMOVED IN THE NEAR FUTURE

auto BipModelPHREEQC(Strings const& substances, BipModelParamsPhreeqc const& params) -> BipModel
{
    return BipModelPhreeqc(substances, params);
}

} // namespace CubicEOS
} // namespace Reaktoro
