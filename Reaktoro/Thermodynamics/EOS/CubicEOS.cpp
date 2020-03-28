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

using AlphaResult = std::tuple<ThermoScalar, ThermoScalar, ThermoScalar>;

/// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT (temperature derivatives) for a given EOS.
auto alpha(CubicEOS::Model type) -> std::function<AlphaResult(const ThermoScalar&, double)>
{
    // The alpha function for van der Waals EOS
    auto alphaVDW = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        ThermoScalar val(1.0);
        ThermoScalar ddt(0.0);
        ThermoScalar d2dt2(0.0);
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Redlich-Kwong EOS
    auto alphaRK = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        ThermoScalar val = 1.0/sqrt(Tr);
        ThermoScalar ddt = -0.5/Tr * val;
        ThermoScalar d2dt2 = -0.5/Tr * (ddt - val/Tr);
        ddt *= Tr.ddT;
        d2dt2 *= Tr.ddT*Tr.ddT;
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Soave-Redlich-Kwong EOS
    auto alphaSRK = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        double m = 0.480 + 1.574*omega - 0.176*omega*omega;
        ThermoScalar sqrtTr = sqrt(Tr);
        ThermoScalar aux_val = 1.0 + m*(1.0 - sqrtTr);
        ThermoScalar aux_ddt = -0.5*m/sqrtTr;
        ThermoScalar aux_d2dt2 = 0.25*m/(Tr*sqrtTr);
        ThermoScalar val = aux_val*aux_val;
        ThermoScalar ddt = 2.0*aux_val*aux_ddt;
        ThermoScalar d2dt2 = 2.0*(aux_ddt*aux_ddt + aux_val*aux_d2dt2);
        ddt *= Tr.ddT;
        d2dt2 *= Tr.ddT*Tr.ddT;
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Peng-Robinson (1978) EOS
    auto alphaPR = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        // Jaubert, J.-N., Vitu, S., Mutelet, F. and Corriou, J.-P., 2005.
        // Extension of the PPR78 model (predictive 1978, Peng–Robinson EOS
        // with temperature dependent kij calculated through a group
        // contribution method) to systems containing aromatic compounds.
        // Fluid Phase Equilibria, 237(1-2), pp.193–211.
        double m = omega <= 0.491 ?
            0.374640 + 1.54226*omega - 0.269920*omega*omega :
            0.379642 + 1.48503*omega - 0.164423*omega*omega + 0.016666*omega*omega*omega;
        ThermoScalar sqrtTr = sqrt(Tr);
        ThermoScalar aux_val = 1.0 + m*(1.0 - sqrtTr);
        ThermoScalar aux_ddt = -0.5*m/sqrtTr;
        ThermoScalar aux_d2dt2 = 0.25*m/(Tr*sqrtTr);
        ThermoScalar val = aux_val*aux_val;
        ThermoScalar ddt = 2.0*aux_val*aux_ddt;
        ThermoScalar d2dt2 = 2.0*(aux_ddt*aux_ddt + aux_val*aux_d2dt2);
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
        ChemicalScalar sca(nspecies);
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

    auto operator()(const ThermoScalar& T, const ThermoScalar& P, const VectorXd& x) -> Result
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
        ThermoVector a(nspecies);
        ThermoVector aT(nspecies);
        ThermoVector aTT(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tc = critical_temperatures[i];
            const double Pc = critical_pressures[i];
            const double omega = acentric_factors[i];
            const double factor = Psi*R*R*(Tc*Tc)/Pc;
            const ThermoScalar Tr = T/Tc;
            ThermoScalar alpha_val, alpha_ddt, alpha_d2dt2;
            std::tie(alpha_val, alpha_ddt, alpha_d2dt2) = alpha(Tr, omega);
            a[i] = factor * alpha_val;
            aT[i] = factor * alpha_ddt;
            aTT[i] = factor * alpha_d2dt2;
        };

        // Calculate the parameters `b` of the cubic equation of state for each species
        Vector b(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tci = critical_temperatures[i];
            const double Pci = critical_pressures[i];
            b[i] = Omega*R*Tci/Pci;
        }

        // Calculate the table of binary interaction parameters and its temperature derivatives
        InteractionParamsResult kres;

        if(calculate_interaction_params)
            kres = calculate_interaction_params(T.val);
        // Calculate the parameter `amix` of the phase and the partial molar parameters `abar` of each species
        ChemicalScalar amix(nspecies);
        ChemicalScalar amixT(nspecies);
        ChemicalScalar amixTT(nspecies);
        VectorXd abar(nspecies);
        VectorXd abarT(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            for(unsigned j = 0; j < nspecies; ++j)
            {
                const ThermoScalar r = kres.k.size() ? ThermoScalar(1.0 - kres.k(i, j)) : ThermoScalar(1.0);
                const ThermoScalar rT = kres.kT.size() ? ThermoScalar(-kres.kT(i, j)) : ThermoScalar(0.0);
                const ThermoScalar rTT = kres.kTT.size() ? ThermoScalar(-kres.kTT(i, j)) : ThermoScalar(0.0);

                const ThermoScalar s = sqrt(a[i]*a[j]);
                const ThermoScalar sT = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
                const ThermoScalar sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;

                const ThermoScalar aij = r*s;
                const ThermoScalar aijT = rT*s + r*sT;
                const ThermoScalar aijTT = rTT*s + 2.0*rT*sT + r*sTT;

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
        ChemicalScalar bmix(nspecies);
        Vector bbar(nspecies);
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
        const ChemicalScalar beta = P*bmix/(R*T);
        const ChemicalScalar betaT = beta * (bmixT/bmix - 1.0/T);

        const ChemicalScalar q = amix/(bmix*R*T);
        const ChemicalScalar qT = q*(amixT/amix - 1.0/T);
        const ChemicalScalar qTT = qT*qT/q + q*(1.0/(T*T) + amixTT/amix - amixT*amixT/(amix*amix));

        // Calculate the coefficients A, B, C of the cubic equation of state
        const ChemicalScalar A = (epsilon + sigma - 1)*beta - 1;
        const ChemicalScalar B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon + sigma - q)*beta;
        const ChemicalScalar C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

        // Calculate the partial temperature derivative of the coefficients A, B, C
        const ChemicalScalar AT = (epsilon + sigma - 1)*betaT;
        const ChemicalScalar BT = 2*(epsilon*sigma - epsilon - sigma)*beta*betaT + qT*beta - (epsilon + sigma - q)*betaT;
        const ChemicalScalar CT = -3*epsilon*sigma*beta*beta*betaT - qT*beta*beta - 2*(epsilon*sigma + q)*beta*betaT;

        // Calculate cubicEOS roots
        const auto roots = cardano(1, A.val, B.val, C.val);

        //verify if all cubicEOS_roots are greater than bmix
        std::vector<double> cubicEOS_roots;
        cubicEOS_roots.reserve(3);
        if (std::get<0>(roots).imag() == 0 && std::get<0>(roots).real() >= bmix.val)
            cubicEOS_roots.push_back(std::get<0>(roots).real());

        if (std::get<1>(roots).imag() == 0 && std::get<1>(roots).real() >= bmix.val)
            cubicEOS_roots.push_back(std::get<1>(roots).real());

        if (std::get<2>(roots).imag() == 0 && std::get<2>(roots).real() >= bmix.val)
            cubicEOS_roots.push_back(std::get<2>(roots).real());

        // All possible Compressibility factor
        const auto cubic_size = cubicEOS_roots.size();

        std::vector<ChemicalScalar> Zs;
        if (cubic_size == 1 || cubic_size == 2)
        {
            //even if cubicEOS_roots has 2 roots, assume that the smallest does not have physical meaning
            Zs.push_back(ChemicalScalar(nspecies, cubicEOS_roots[0]));
        }
        else
        {
            if (cubic_size != 3) {
                Exception exception;
                exception.error << "Could not calculate the cubic equation of state.";
                exception.reason << "Logic error: it was expected Z roots of size 3, but got: " << Zs.size();
                RaiseError(exception);
            }
            Zs.push_back(ChemicalScalar(nspecies, cubicEOS_roots[0]));  // Z_max
            Zs.push_back(ChemicalScalar(nspecies, cubicEOS_roots[2]));  // Z_min
        }

        // Selecting compressibility factor - Z_liq < Z_gas
        ChemicalScalar Z(nspecies);
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

        ChemicalScalar& V = result.molar_volume;
        ChemicalScalar& G_res = result.residual_molar_gibbs_energy;
        ChemicalScalar& H_res = result.residual_molar_enthalpy;
        ChemicalScalar& Cp_res = result.residual_molar_heat_capacity_cp;
        ChemicalScalar& Cv_res = result.residual_molar_heat_capacity_cv;
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
        const ChemicalScalar ZT = -(AT*Z*Z + BT*Z + CT)/(3*Z*Z + 2*A*Z + B);

        // Calculate the integration factor I and its temperature derivative IT
        ChemicalScalar I;
        if(epsilon != sigma) I = log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon);
                        else I = beta/(Z + epsilon*beta);

        // Calculate the temperature derivative IT of the integration factor I
        ChemicalScalar IT;
        if(epsilon != sigma) IT = ((ZT + sigma*betaT)/(Z + sigma*beta) - (ZT + epsilon*betaT)/(Z + epsilon*beta))/(sigma - epsilon);
                        else IT = I*(betaT/beta - (ZT + epsilon*betaT)/(Z + epsilon*beta));


        // Calculate the partial molar Zi for each species
        V = Z*R*T/P;
        G_res = R*T*(Z - 1 - log(Z - beta) - q*I);
        H_res = R*T*(Z - 1 + T*qT*I);
        Cp_res = R*T*(ZT + qT*I + T*qTT + T*qT*IT) + H_res/T;

        const ChemicalScalar dPdT = P*(1.0/T + ZT/Z);
        const ChemicalScalar dVdT = V*(1.0/T + ZT/Z);

        Cv_res = Cp_res - T * dPdT*dVdT + R;


        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double bi = bbar[i];
            const ThermoScalar betai = P*bi/(R*T);
            const ChemicalScalar ai = abar[i];
            const ChemicalScalar aiT = abarT[i];
            const ChemicalScalar qi = q*(1 + ai/amix - bi/bmix);
            const ChemicalScalar qiT = qi*qT/q + q*(aiT - ai*amixT/amix)/amix;
            const ThermoScalar Ai = (epsilon + sigma - 1.0)*betai - 1.0;
            const ChemicalScalar Bi = (epsilon*sigma - epsilon - sigma)*(2*beta*betai - beta*beta) - (epsilon + sigma - q)*(betai - beta) - (epsilon + sigma - qi)*beta;
            const ChemicalScalar Ci = -3*sigma*epsilon*beta*beta*betai + 2*epsilon*sigma*beta*beta*beta - (epsilon*sigma + qi)*beta*beta - 2*(epsilon*sigma + q)*(beta*betai - beta*beta);
            const ChemicalScalar Zi = -(Ai*Z*Z + (Bi + B)*Z + Ci + 2*C)/(3*Z*Z + 2*A*Z + B);

            ChemicalScalar Ii;
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

/// Sanity check free function to verify if BIPs matrices have proper dimensions. Considering that the phase has
/// n species, the BIP matricies k, kT and kTT should have (n, n) as dimensions.
/// @see CubicEOS::setInteractionParamsFunction
auto sanityCheckInteractionParamsFunction(const unsigned& nspecies, const CubicEOS::InteractionParamsFunction& func) -> void
{
    auto T_for_sanity_check = 273.0;
    auto bips = func(T_for_sanity_check);

    // Check k's dimensions
    auto k_size = bips.k.size();
    if (k_size > 0) {
        auto k_num_of_rows = bips.k.rows();
        auto k_num_of_cols = bips.k.cols();
        Assert(k_size == nspecies * nspecies && k_num_of_cols == nspecies && k_num_of_rows == nspecies,
            "Could not set the binary interaction parameters (k) in the CubicEOS.",
            "Dimension mismatch between given BIP matrix and number of species.");

        // Check k's symmetry
        for (unsigned i = 0; i < k_num_of_rows; i++){
            for (auto j = i + 1; j < k_num_of_cols; j++){
                Assert(bips.k(i, j) == bips.k(j, i),
                    "BIPs matrix k is not symmetric.", "Check your k BIPs matrix input."
                );
            }
        }
    }
    else {
        Exception exception;
        exception.error << "Invalid function to calculate the BIPs matrix k.";
        exception.reason << "Logic error: the function is unable to return a BIP matrix k.";
        RaiseError(exception);
    }

    // Check kT's dimensions
    auto kT_size = bips.kT.size();
    if (kT_size > 0) {  // if kT is provided
        auto kT_num_of_rows = bips.kT.rows();
        auto kT_num_of_cols = bips.kT.cols();
        Assert(kT_size == nspecies * nspecies && kT_num_of_cols == nspecies && kT_num_of_rows == nspecies,
            "Could not set the binary interaction parameters (kT) in the CubicEOS.",
            "Dimension mismatch between given BIP matrix and number of species.");

        // Check kT's symmetry
        for (unsigned i = 0; i < kT_num_of_rows; i++){
            for (auto j = i + 1; j < kT_num_of_cols; j++){
                Assert(bips.kT(i, j) == bips.kT(j, i),
                    "BIPs matrix kT is not symmetric.", "Check your kT BIPs matrix input.");
            }
        }
    }

    // Check kTT's dimensions
    auto kTT_size = bips.kTT.size();
    if (kTT_size > 0) {  // if kTT is provided
        auto kTT_num_of_rows = bips.kTT.rows();
        auto kTT_num_of_cols = bips.kTT.cols();
        Assert(kTT_size == nspecies * nspecies && kTT_num_of_cols == nspecies && kTT_num_of_rows == nspecies,
            "Could not set the binary interaction parameters (kTT) in the CubicEOS.",
            "Dimension mismatch between given BIP matrix and number of species.");

        // Check kTT's symmetry
        for (unsigned i = 0; i < kTT_num_of_rows; i++){
            for (auto j = i + 1; j < kTT_num_of_cols; j++){
                Assert(bips.kT(i, j) == bips.kT(j, i),
                    "BIPs matrix kT is not symmetric.", "Check your kT BIPs matrix input."
                );
            }
        }
    }
}

CubicEOS::CubicEOS(unsigned nspecies, CubicEOS::Params params)
: pimpl(new Impl(nspecies))
{
    pimpl->model = params.model;
    pimpl->phase_identification_method = params.phase_identification_method;
    if(params.binary_interaction_values) {
        sanityCheckInteractionParamsFunction(nspecies, params.binary_interaction_values);
        setInteractionParamsFunction(params.binary_interaction_values);
    }
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

auto CubicEOS::operator()(const ThermoScalar& T, const ThermoScalar& P, const VectorXd& x) -> Result
{
    return pimpl->operator()(T, P, x);
}

} // namespace Reaktoro
