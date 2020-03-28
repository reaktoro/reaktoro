// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Math/Roots.hpp>
#include <Reaktoro/Thermodynamics/EOS/PhaseIdentification.hpp>

namespace Reaktoro {
namespace internal {

using AlphaResult = std::tuple<real, real, real>;

/// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT (temperature derivatives) for a given EOS.
auto alpha(CubicEOS::Model type) -> std::function<AlphaResult(const real&, double)>
{
    // The alpha function for van der Waals EOS
    auto alphaVDW = [](const real& Tr, double omega) -> AlphaResult
    {
        real val(1.0);
        real ddt(0.0);
        real d2dt2(0.0);
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Redlich-Kwong EOS
    auto alphaRK = [](const real& Tr, double omega) -> AlphaResult
    {
        real val = 1.0/sqrt(Tr);
        real ddt = -0.5/Tr * val;
        real d2dt2 = -0.5/Tr * (ddt - val/Tr);
        ddt *= Tr.ddT;
        d2dt2 *= Tr.ddT*Tr.ddT;
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Soave-Redlich-Kwong EOS
    auto alphaSRK = [](const real& Tr, double omega) -> AlphaResult
    {
        double m = 0.480 + 1.574*omega - 0.176*omega*omega;
        real sqrtTr = sqrt(Tr);
        real aux_val = 1.0 + m*(1.0 - sqrtTr);
        real aux_ddt = -0.5*m/sqrtTr;
        real aux_d2dt2 = 0.25*m/(Tr*sqrtTr);
        real val = aux_val*aux_val;
        real ddt = 2.0*aux_val*aux_ddt;
        real d2dt2 = 2.0*(aux_ddt*aux_ddt + aux_val*aux_d2dt2);
        ddt *= Tr.ddT;
        d2dt2 *= Tr.ddT*Tr.ddT;
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Peng-Robinson (1978) EOS
    auto alphaPR = [](const real& Tr, double omega) -> AlphaResult
    {
        // Jaubert, J.-N., Vitu, S., Mutelet, F. and Corriou, J.-P., 2005.
        // Extension of the PPR78 model (predictive 1978, Peng–Robinson EOS
        // with temperature dependent kij calculated through a group
        // contribution method) to systems containing aromatic compounds.
        // Fluid Phase Equilibria, 237(1-2), pp.193–211.
        double m = omega <= 0.491 ?
            0.374640 + 1.54226*omega - 0.269920*omega*omega :
            0.379642 + 1.48503*omega - 0.164423*omega*omega + 0.016666*omega*omega*omega;
        real sqrtTr = sqrt(Tr);
        real aux_val = 1.0 + m*(1.0 - sqrtTr);
        real aux_ddt = -0.5*m/sqrtTr;
        real aux_d2dt2 = 0.25*m/(Tr*sqrtTr);
        real val = aux_val*aux_val;
        real ddt = 2.0*aux_val*aux_ddt;
        real d2dt2 = 2.0*(aux_ddt*aux_ddt + aux_val*aux_d2dt2);
        ddt *= Tr.ddT;
        d2dt2 *= Tr.ddT*Tr.ddT;
        return std::make_tuple(val, ddt, d2dt2);
    };

    switch(type)
    {
        case CubicEOS::VanDerWaals: return alphaVDW;
        case CubicEOS::RedlichKwong: return alphaRK;
        case CubicEOS::SoaveRedlichKwong: return alphaSRK;
        case CubicEOS::PengRobinson: return alphaPR;
        default: return alphaPR;
    }
}

auto sigma(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 0.0;
        case CubicEOS::RedlichKwong: return 1.0;
        case CubicEOS::SoaveRedlichKwong: return 1.0;
        case CubicEOS::PengRobinson: return 1.0 + 1.4142135623730951;
        default: return 1.0 + 1.4142135623730951;
    }
}

auto epsilon(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 0.0;
        case CubicEOS::RedlichKwong: return 0.0;
        case CubicEOS::SoaveRedlichKwong: return 0.0;
        case CubicEOS::PengRobinson: return 1.0 - 1.4142135623730951;
        default: return 1.0 - 1.4142135623730951;
    }
}

auto Omega(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 1.0/8.0;
        case CubicEOS::RedlichKwong: return 0.08664;
        case CubicEOS::SoaveRedlichKwong: return 0.08664;
        case CubicEOS::PengRobinson: return 0.0777960739;
        default: return 0.0777960739;
    }
}

auto Psi(CubicEOS::Model type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 27.0/64.0;
        case CubicEOS::RedlichKwong: return 0.42748;
        case CubicEOS::SoaveRedlichKwong: return 0.42748;
        case CubicEOS::PengRobinson: return 0.457235529;
        default: return 0.457235529;
    }
}

} // namespace internal

struct CubicEOS::Impl
{
    /// The number of species in the phase.
    unsigned nspecies;

    /// The flag that indicates if the phase is vapor (false means liquid instead).
    bool isvapor = true;

    /// The type of the cubic equation of state.
    CubicEOS::Model model = CubicEOS::PengRobinson;

    /// The type of phase identification method that is going to be used
    PhaseIdentificationMethod phase_identification_method = PhaseIdentificationMethod::None;

    /// The critical temperatures of the species (in units of K).
    std::vector<double> critical_temperatures;

    /// The critical pressures of the species (in units of Pa).
    std::vector<double> critical_pressures;

    /// The acentric factors of the species.
    std::vector<double> acentric_factors;

    /// The function that calculates the interaction parameters kij and its temperature derivatives.
    InteractionParamsFunction calculate_interaction_params;

    /// The result with thermodynamic properties calculated from the cubic equation of state
    Result result;

    /// Construct a CubicEOS::Impl instance.
    Impl(unsigned nspecies)
    : nspecies(nspecies)
    {
        // Initialize the dimension of the chemical scalar quantities
        real sca(nspecies);
        result.molar_volume = sca;
        result.residual_molar_gibbs_energy = sca;
        result.residual_molar_enthalpy = sca;
        result.residual_molar_heat_capacity_cp = sca;
        result.residual_molar_heat_capacity_cv = sca;

        // Initialize the dimension of the chemical vector quantities
        VectorXd vec(nspecies);
        result.partial_molar_volumes = vec;
        result.residual_partial_molar_enthalpies = vec;
        result.residual_partial_molar_gibbs_energies = vec;
        result.ln_fugacity_coefficients = vec;
    }

    auto operator()(const real& T, const real& P, const VectorXd& x) -> Result
    {
        // Check if the mole fractions are zero or non-initialized
        if(x.val.size() == 0 || min(x.val) <= 0.0)
            return Result(nspecies); // result with zero values

        // Auxiliary variables
        const double R = universalGasConstant;
        const double Psi = internal::Psi(model);
        const double Omega = internal::Omega(model);
        const double epsilon = internal::epsilon(model);
        const double sigma = internal::sigma(model);
        const auto alpha = internal::alpha(model);

        // Calculate the parameters `a` of the cubic equation of state for each species
        VectorXr a(nspecies);
        VectorXr aT(nspecies);
        VectorXr aTT(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tc = critical_temperatures[i];
            const double Pc = critical_pressures[i];
            const double omega = acentric_factors[i];
            const double factor = Psi*R*R*(Tc*Tc)/Pc;
            const real Tr = T/Tc;
            real alpha_val, alpha_ddt, alpha_d2dt2;
            std::tie(alpha_val, alpha_ddt, alpha_d2dt2) = alpha(Tr, omega);
            a[i] = factor * alpha_val;
            aT[i] = factor * alpha_ddt;
            aTT[i] = factor * alpha_d2dt2;
        };

        // Calculate the parameters `b` of the cubic equation of state for each species
        VectorXr b(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tci = critical_temperatures[i];
            const double Pci = critical_pressures[i];
            b[i] = Omega*R*Tci/Pci;
        }

        // Calculate the table of binary interaction parameters and its temperature derivatives
        InteractionParamsResult kres;
        InteractionParamsArgs kargs{T, a, aT, aTT, b};

        if(calculate_interaction_params)
            kres = calculate_interaction_params(kargs);
        // Calculate the parameter `amix` of the phase and the partial molar parameters `abar` of each species
        real amix(nspecies);
        real amixT(nspecies);
        real amixTT(nspecies);
        VectorXd abar(nspecies);
        VectorXd abarT(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            for(unsigned j = 0; j < nspecies; ++j)
            {
                const real r = kres.k.empty() ? real(1.0) : 1.0 - kres.k[i][j];
                const real rT = kres.kT.empty() ? real(0.0) : -kres.kT[i][j];
                const real rTT = kres.kTT.empty() ? real(0.0) : -kres.kTT[i][j];

                const real s = sqrt(a[i]*a[j]);
                const real sT = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
                const real sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;

                const real aij = r*s;
                const real aijT = rT*s + r*sT;
                const real aijTT = rTT*s + 2.0*rT*sT + r*sTT;

                amix += x[i] * x[j] * aij;
                amixT += x[i] * x[j] * aijT;
                amixTT += x[i] * x[j] * aijTT;

                abar[i] += 2 * x[j] * aij;
                abarT[i] += 2 * x[j] * aijT;
            }
        }

        // Finalize the calculation of `abar` and `abarT`
        for(unsigned i = 0; i < nspecies; ++i)
        {
            abar[i] -= amix;
            abarT[i] -= amixT;
        }

        // Calculate the parameter `bmix` of the cubic equation of state
        real bmix(nspecies);
        VectorXr bbar(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tci = critical_temperatures[i];
            const double Pci = critical_pressures[i];
            bbar[i] = Omega*R*Tci/Pci;
            bmix += x[i] * bbar[i];
        }

        // Calculate the temperature derivative of `bmix`
        const double bmixT = 0.0; // no temperature dependence

        // Calculate auxiliary quantities `beta` and `q`
        const real beta = P*bmix/(R*T);
        const real betaT = beta * (bmixT/bmix - 1.0/T);

        const real q = amix/(bmix*R*T);
        const real qT = q*(amixT/amix - 1.0/T);
        const real qTT = qT*qT/q + q*(1.0/(T*T) + amixTT/amix - amixT*amixT/(amix*amix));

        // Calculate the coefficients A, B, C of the cubic equation of state
        const real A = (epsilon + sigma - 1)*beta - 1;
        const real B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon + sigma - q)*beta;
        const real C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

        // Calculate the partial temperature derivative of the coefficients A, B, C
        const real AT = (epsilon + sigma - 1)*betaT;
        const real BT = 2*(epsilon*sigma - epsilon - sigma)*beta*betaT + qT*beta - (epsilon + sigma - q)*betaT;
        const real CT = -3*epsilon*sigma*beta*beta*betaT - qT*beta*beta - 2*(epsilon*sigma + q)*beta*betaT;

        // Calculate cubicEOS roots using cardano's method
        auto cubicEOS_roots = realRoots(cardano(1, A.val, B.val, C.val));

        // All possible Compressibility factor
        std::vector<real> Zs;
        if (cubicEOS_roots.size() == 1)
        {
            Zs.push_back(real(nspecies, cubicEOS_roots[0]));
        }
        else
        {
            if (cubicEOS_roots.size() != 3) {
                Exception exception;
                exception.error << "Could not calculate the cubic equation of state.";
                exception.reason << "Logic error: it was expected Z roots of size 3, but got: " << Zs.size();
                RaiseError(exception);
            }
            Zs.push_back(real(nspecies, cubicEOS_roots[0]));  // Z_max
            Zs.push_back(real(nspecies, cubicEOS_roots[2]));  // Z_min
        }

        // Selecting compressibility factor - Z_liq < Z_gas
        real Z(nspecies);
        if (isvapor)
            Z.val = *std::max_element(cubicEOS_roots.begin(), cubicEOS_roots.end());
        else
            Z.val = *std::min_element(cubicEOS_roots.begin(), cubicEOS_roots.end());

        auto input_phase_type = isvapor ? PhaseType::Gas : PhaseType::Liquid;
        auto identified_phase_type = input_phase_type;

        switch (phase_identification_method)
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

        if (identified_phase_type != input_phase_type)
        {
            // Since the phase is identified as different than the expect input phase type, it is
            // deemed inappropriate. Artificially high values are configured for fugacities, so that
            // this condition is "removed" by the optimizer.
            result.molar_volume = 0.0;
            result.residual_molar_gibbs_energy = 0.0;
            result.residual_molar_enthalpy = 0.0;
            result.residual_molar_heat_capacity_cp = 0.0;
            result.residual_molar_heat_capacity_cv = 0.0;
            result.partial_molar_volumes.fill(0.0);
            result.residual_partial_molar_gibbs_energies.fill(0.0);
            result.residual_partial_molar_enthalpies.fill(0.0);
            result.ln_fugacity_coefficients.fill(100.0);
            return result;
        }

        real& V = result.molar_volume;
        real& G_res = result.residual_molar_gibbs_energy;
        real& H_res = result.residual_molar_enthalpy;
        real& Cp_res = result.residual_molar_heat_capacity_cp;
        real& Cv_res = result.residual_molar_heat_capacity_cv;
        VectorXd& Vi = result.partial_molar_volumes;
        VectorXd& Gi_res = result.residual_partial_molar_gibbs_energies;
        VectorXd& Hi_res = result.residual_partial_molar_enthalpies;
        VectorXd& ln_phi = result.ln_fugacity_coefficients;

        // Calculate the partial derivatives of Z (dZdT, dZdP, dZdn)
        const double factor = -1.0/(3*Z.val*Z.val + 2*A.val*Z.val + B.val);
        Z.ddT = factor * (A.ddT*Z.val*Z.val + B.ddT*Z.val + C.ddT);
        Z.ddP = factor * (A.ddP*Z.val*Z.val + B.ddP*Z.val + C.ddP);
        for(unsigned i = 0; i < nspecies; ++i)
            Z.ddn[i] = factor * (A.ddn[i]*Z.val*Z.val + B.ddn[i]*Z.val + C.ddn[i]);

        // Calculate the partial temperature derivative of Z
        const real ZT = -(AT*Z*Z + BT*Z + CT)/(3*Z*Z + 2*A*Z + B);

        // Calculate the integration factor I and its temperature derivative IT
        real I;
        if(epsilon != sigma) I = log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon);
                        else I = beta/(Z + epsilon*beta);

        // Calculate the temperature derivative IT of the integration factor I
        real IT;
        if(epsilon != sigma) IT = ((ZT + sigma*betaT)/(Z + sigma*beta) - (ZT + epsilon*betaT)/(Z + epsilon*beta))/(sigma - epsilon);
                        else IT = I*(betaT/beta - (ZT + epsilon*betaT)/(Z + epsilon*beta));


        // Calculate the partial molar Zi for each species
        V = Z*R*T/P;
        G_res = R*T*(Z - 1 - log(Z - beta) - q*I);
        H_res = R*T*(Z - 1 + T*qT*I);
        Cp_res = R*T*(ZT + qT*I + T*qTT + T*qT*IT) + H_res/T;

        const real dPdT = P*(1.0/T + ZT/Z);
        const real dVdT = V*(1.0/T + ZT/Z);

        Cv_res = Cp_res - T * dPdT*dVdT + R;


        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double bi = bbar[i];
            const real betai = P*bi/(R*T);
            const real ai = abar[i];
            const real aiT = abarT[i];
            const real qi = q*(1 + ai/amix - bi/bmix);
            const real qiT = qi*qT/q + q*(aiT - ai*amixT/amix)/amix;
            const real Ai = (epsilon + sigma - 1.0)*betai - 1.0;
            const real Bi = (epsilon*sigma - epsilon - sigma)*(2*beta*betai - beta*beta) - (epsilon + sigma - q)*(betai - beta) - (epsilon + sigma - qi)*beta;
            const real Ci = -3*sigma*epsilon*beta*beta*betai + 2*epsilon*sigma*beta*beta*beta - (epsilon*sigma + qi)*beta*beta - 2*(epsilon*sigma + q)*(beta*betai - beta*beta);
            const real Zi = -(Ai*Z*Z + (Bi + B)*Z + Ci + 2*C)/(3*Z*Z + 2*A*Z + B);

            real Ii;
            if(epsilon != sigma) Ii = I + ((Zi + sigma*betai)/(Z + sigma*beta) - (Zi + epsilon*betai)/(Z + epsilon*beta))/(sigma - epsilon);
                            else Ii = I * (1 + betai/beta - (Zi + epsilon*betai)/(Z + epsilon*beta));

            Vi[i] = R*T*Zi/P;
            Gi_res[i] = R*T*(Zi - (Zi - betai)/(Z - beta) - log(Z - beta) - qi*I - q*Ii + q*I);
            Hi_res[i] = R*T*(Zi - 1 + T*(qiT*I + qT*Ii - qT*I));
            ln_phi[i] = Gi_res[i]/(R*T);
        }

        return result;
    }
};

CubicEOS::Result::Result()
{}

CubicEOS::Result::Result(unsigned nspecies)
: molar_volume(nspecies),
  residual_molar_gibbs_energy(nspecies),
  residual_molar_enthalpy(nspecies),
  residual_molar_heat_capacity_cp(nspecies),
  residual_molar_heat_capacity_cv(nspecies),
  partial_molar_volumes(nspecies),
  residual_partial_molar_gibbs_energies(nspecies),
  residual_partial_molar_enthalpies(nspecies),
  ln_fugacity_coefficients(nspecies)
{}

CubicEOS::CubicEOS(unsigned nspecies, CubicEOS::Params params)
: pimpl(new Impl(nspecies))
{
    pimpl->model = params.model;
    pimpl->phase_identification_method = params.phase_identification_method;
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

auto CubicEOS::setModel(Model model) -> void
{
    pimpl->model = model;
}

auto CubicEOS::setPhaseAsLiquid() -> void
{
    pimpl->isvapor = false;
}

auto CubicEOS::setPhaseAsVapor() -> void
{
    pimpl->isvapor = true;
}

auto CubicEOS::setCriticalTemperatures(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the critical "
        "temperatures of the species in the CubicEOS object.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " were given.");

    Assert(minValue(values) > 0, "Cannot set the critical temperatures of the species "
        "in the CubicEOS object.", "Expecting non-zero critical "
        "temperatures of the gases.");

    pimpl->critical_temperatures = values;
}

auto CubicEOS::setCriticalPressures(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the critical "
        "pressures of the species in the CubicEOS object.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " were given.");

    Assert(minValue(values) > 0, "Cannot set the critical pressures of the species "
        "in the CubicEOS object.", "Expecting non-zero critical "
        "pressures of the gases.");

    pimpl->critical_pressures = values;
}

auto CubicEOS::setAcentricFactors(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the acentric "
        "factors of the species in CubicEOS.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " values were given.");

    pimpl->acentric_factors = values;
}

auto CubicEOS::setInteractionParamsFunction(const InteractionParamsFunction& func) -> void
{
    pimpl->calculate_interaction_params = func;
}

auto CubicEOS::operator()(const real& T, const real& P, const VectorXd& x) -> Result
{
    return pimpl->operator()(T, P, x);
}

} // namespace Reaktoro
