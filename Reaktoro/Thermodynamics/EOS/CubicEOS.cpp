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

#include "CubicEOS.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/TableUtils.hpp>
#include <Reaktoro/Math/Roots.hpp>
#include <Reaktoro/Thermodynamics/EOS/PhaseIdentification.hpp>

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

using std::log;
using std::sqrt;

namespace internal {

using AlphaResult = std::tuple<real, real, real>;

/// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT for a given EOS.
auto alpha(CubicEOSModel type) -> std::function<AlphaResult(real, real, real)>
{
    // The alpha function for van der Waals EOS (see Table 3.1 of Smith et al. 2017)
    auto alphaVDW = [](real T, real Tc, real omega) -> AlphaResult
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
        const real alphaT = -0.5/Tr*alpha * TrT; // === -0.5/(Tr*sqrt(Tr)) * TrT === -0.5/Tr*alpha * TrT
        const real alphaTT = 0.5/Tr*(alpha/Tr * TrT - alphaT) * TrT; // === (0.5/Tr*alpha/Tr * TrT - 0.5/Tr*alphaT) * TrT === 0.5/Tr*(alpha/Tr * TrT - alphaT) * TrT
        return { alpha, alphaT, alphaTT };
    };

    // The alpha function for Soave-Redlich-Kwong EOS
    auto alphaSRK = [](real Tr, real TrT, real omega) -> AlphaResult
    {
        const real m = 0.480 + 1.574*omega - 0.176*omega*omega;
        const real sqrtTr = sqrt(Tr);
        const real sqrtTrT = 0.5/sqrtTr * TrT; // === 0.5/sqrt(Tr) * TrT
        const real sqrtTrTT = -0.5/Tr * sqrtTrT * TrT; // === -0.5/(sqrtTr*sqrtTr) * sqrtTrT * TrT
        const real aux = 1.0 + m*(1.0 - sqrtTr);
        const real auxT = -m*sqrtTrT;
        const real auxTT = -m*sqrtTrTT;
        const real alpha = aux*aux;
        const real alphaT = 2.0*aux*auxT;
        const real alphaTT = 2.0*(auxT*auxT + aux*auxTT);
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
        const real sqrtTrT = 0.5/sqrtTr * TrT; // === 0.5/sqrt(Tr) * TrT
        const real sqrtTrTT = -0.5/Tr * sqrtTrT * TrT; // === -0.5/(sqrtTr*sqrtTr) * sqrtTrT * TrT
        const real aux = 1.0 + m*(1.0 - sqrtTr);
        const real auxT = -m*sqrtTrT;
        const real auxTT = -m*sqrtTrTT;
        const real alpha = aux*aux;
        const real alphaT = 2.0*aux*auxT;
        const real alphaTT = 2.0*(auxT*auxT + aux*auxTT);
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

} // namespace internal

struct CubicEOS::Impl
{
    /// The number of species in the phase.
    unsigned nspecies;

    /// The fluid type for which the equation of state should be confifured.
    CubicEOSFluidType fluidtype = CubicEOSFluidType::Vapor;

    /// The type of the cubic equation of state.
    CubicEOSModel model = CubicEOSModel::PengRobinson;

    /// The type of phase identification method that is going to be used
    PhaseIdentificationMethod phase_identification_method = PhaseIdentificationMethod::None;

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
      fluidtype(args.fluidtype),
      model(args.model),
      phase_identification_method(args.phase_identification_method),
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
        const auto R = universalGasConstant;
        const auto Psi = internal::Psi(model);
        const auto Omega = internal::Omega(model);
        const auto epsilon = internal::epsilon(model);
        const auto sigma = internal::sigma(model);
        const auto alphafn = internal::alpha(model);

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
        };

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
                const real r   = ip.k.size()   ? 1.0 - ip.k(i, j)   : 1.0;
                const real rT  = ip.kT.size()  ?     - ip.kT(i, j)  : 0.0;
                const real rTT = ip.kTT.size() ?     - ip.kTT(i, j) : 0.0;

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

        // Calculate cubicEOS roots using cardano's method
        auto cubicEOS_roots = realRoots(cardano(A, B, C));

        // All possible Compressibility factor
        std::vector<real> Zs;
        if (cubicEOS_roots.size() == 1)
        {
            Zs.push_back(cubicEOS_roots[0]);
        }
        else
        {
            if (cubicEOS_roots.size() != 3) {
                Exception exception;
                exception.error << "Could not calculate the cubic equation of state.";
                exception.reason << "Logic error: it was expected Z roots of size 3, but got: " << Zs.size();
                RaiseError(exception);
            }
            Zs.push_back(cubicEOS_roots[0]);  // Z_max
            Zs.push_back(cubicEOS_roots[2]);  // Z_min
        }

        // Selecting compressibility factor - Z_liq < Z_gas
        real Z = {};
        if (fluidtype == CubicEOSFluidType::Vapor)
            Z = *std::max_element(cubicEOS_roots.begin(), cubicEOS_roots.end());
        else
            Z = *std::min_element(cubicEOS_roots.begin(), cubicEOS_roots.end());

        auto identified_phase_type = fluidtype;

        switch(phase_identification_method)
        {
        case PhaseIdentificationMethod::None:
            // `identified_phase_type` is already `input_phase_type`, keep it this way
            break;

        case PhaseIdentificationMethod::VolumeMethod:
            identified_phase_type = identifyPhaseUsingVolume(T, P, Z, bmix);
            break;

        case PhaseIdentificationMethod::IsothermalCompressibilityMethods:
            identified_phase_type = identifyPhaseUsingIsothermalCompressibility(T, P, Z);
            break;

        case PhaseIdentificationMethod::GibbsEnergyAndEquationOfStateMethod:
            identified_phase_type = identifyPhaseUsingGibbsEnergyAndEos(
                P, T, amix, bmix, A, B, C, Zs, epsilon, sigma);
            break;

        default:
            throw std::logic_error("CubicEOS received an unexpected phaseIdentificationMethod");
        }

        if(identified_phase_type != fluidtype)
        {
            // Since the phase is identified as different than the expect input phase type, it is
            // deemed inappropriate. Artificially high values are configured for fugacities, so that
            // this condition is "removed" by the optimizer.
            props.ln_phi.fill(100.0);
            return;
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

        // Auxiliary references
        auto& [Vres, VresT, VresP, Gres, Hres, Cpres, Cvres, ln_phi] = props;

        //=========================================================================================
        // Calculate the ideal volume properties of the phase
        //=========================================================================================
        const real V0  =  R*T/P;
        const real V0T =  V0/T;
        const real V0P = -V0/P;

        //=========================================================================================
        // Calculate the real volume properties of the phase
        //=========================================================================================
        const real V  = Z*V0;
        const real VT = ZT*V0 + Z*V0T;
        const real VP = ZP*V0 + Z*V0P;

        //=========================================================================================
        // Calculate the residual properties of the phase
        //=========================================================================================
        Vres  = V - Vres;
        VresT = VT - VresT;
        VresP = VP - VresP;
        Gres  = R*T*(Z - 1 - log(Z - beta) - q*I); // from Eq. (13.74) of Smith et al. (2017)
        Hres  = R*T*(Z - 1 + T*qT*I); // equation after Eq. (13.74), but using T*qT instead of Tr*qTr, which is equivalent
        Cpres = Hres/T + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT); // from Eq. (2.19), Cp(res) := (dH(res)/dT)P === R*(Z - 1 + T*qT*I) + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT) = H_res/T + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT)
        Cvres = Cpres + R + T*VT*VT/VP; // from Cv = Cp + T*(VT*VT)/VP and Cv0 = Cp0 - R, Cv(res) = Cp(res) + R + T*(VT*VT)/VP

        //=========================================================================================
        // Calculate the fugacity coefficients for each species
        //=========================================================================================
        for(auto k = 0; k < nspecies; ++k)
        {
            const real betak = bbar[k]/b[k] * beta;
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

auto CubicEOS::numSpecies() const -> unsigned
{
    return pimpl->nspecies;
}

auto CubicEOS::setModel(CubicEOSModel model) -> void
{
    pimpl->model = model;
}

auto CubicEOS::setFluidType(CubicEOSFluidType fluidtype) -> void
{
    pimpl->fluidtype = fluidtype;
}

auto CubicEOS::setInteractionParamsFunction(const CubicEOSInteractionParamsFn& func) -> void
{
    pimpl->interaction_params_fn = func;
}

auto CubicEOS::setStablePhaseIdentificationMethod(const PhaseIdentificationMethod& method) -> void
{
    pimpl->phase_identification_method = method;
}

auto CubicEOS::compute(CubicEOSProps& props, real T, real P, ArrayXrConstRef x) -> void
{
    return pimpl->compute(props, T, P, x);
}

} // namespace Reaktoro
