// Reaktoro is a C++ library for computational reaction modelling.
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

#include "CubicEOS.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/TableUtils.hpp>
#include <Reaktoro/Math/Roots.hpp>

namespace Reaktoro {
namespace internal {

using AlphaResult = std::tuple<ThermoScalar, ThermoScalar, ThermoScalar>;

/// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT (temperature derivatives) for a given EOS.
auto alpha(CubicEOS::Type type) -> std::function<AlphaResult(const ThermoScalar&, double)>
{
    // The alpha function for van der Walls EOS
    auto alphaVDW = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        ThermoScalar val = 1.0;
        ThermoScalar ddt = 0.0;
        ThermoScalar d2dt2 = 0.0;
        return std::make_tuple(val, ddt, d2dt2);
    };

    // The alpha function for Redlich-Kwong EOS
    auto alphaRK = [](const ThermoScalar& Tr, double omega) -> AlphaResult
    {
        ThermoScalar val = 1.0/sqrt(Tr);
        ThermoScalar ddt = -0.5/Tr * val;
        ThermoScalar d2dt2 = -0.5/Tr * (ddt - val/Tr);
        ddt *= Tr.ddt;
        d2dt2 *= Tr.ddt*Tr.ddt;
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
        ThermoScalar ddt = 2*aux_val*aux_ddt;
        ThermoScalar d2dt2 = 2*(aux_ddt*aux_ddt + aux_val*aux_d2dt2);
        ddt *= Tr.ddt;
        d2dt2 *= Tr.ddt*Tr.ddt;
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
        ThermoScalar ddt = 2*aux_val*aux_ddt;
        ThermoScalar d2dt2 = 2*(aux_ddt*aux_ddt + aux_val*aux_d2dt2);
        ddt *= Tr.ddt;
        d2dt2 *= Tr.ddt*Tr.ddt;
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

auto sigma(CubicEOS::Type type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 0.0;
        case CubicEOS::RedlichKwong: return 1.0;
        case CubicEOS::SoaveRedlichKwong: return 1.0;
        case CubicEOS::PengRobinson: return 1.0 + 1.41421356237;
        default: return 1.0 + 1.41421356237;
    }
}

auto epsilon(CubicEOS::Type type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 0.0;
        case CubicEOS::RedlichKwong: return 0.0;
        case CubicEOS::SoaveRedlichKwong: return 0.0;
        case CubicEOS::PengRobinson: return 1.0 - 1.41421356237;
        default: return 1.0 - 1.41421356237;
    }
}

auto Omega(CubicEOS::Type type) -> double
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

auto Psi(CubicEOS::Type type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 27.0/64.0;
        case CubicEOS::RedlichKwong: return 0.42748;
        case CubicEOS::SoaveRedlichKwong: return 0.42748;
        case CubicEOS::PengRobinson: return 0.457235529;
        default: return 0.45724;
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
    CubicEOS::Type eostype = CubicEOS::PengRobinson;

    /// The critical temperatures of the species (in units of K).
    std::vector<double> critical_temperatures;

    /// The critical pressures of the species (in units of Pa).
    std::vector<double> critical_pressures;

    /// The acentric factors of the species.
    std::vector<double> acentric_factors;

    /// The function that calculates the interaction parameters kij and its temperature derivatives.
    InteractionParamsFunction calculate_interaction_params;

//    /// The pairs of species indices for the binary interaction parameters `kij`
//    std::vector<std::tuple<Index, Index>> ij;
//
//    /// The binary intectaction parameters `kij`. The number of interaction parameters
//    /// per pair of species are either **one** or **three**. If three parameters are given,
//    /// these are used in the following correlation equation for temperature:
//    std::vector<std::vector<double>> kij;

    /// Construct a CubicEOS::Impl instance.
    Impl(unsigned nspecies)
    : nspecies(nspecies)
    {}

    auto operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> Result
    {
        // Auxiliary variables
        const double R = universalGasConstant;
        const double Psi = internal::Psi(eostype);
        const double Omega = internal::Omega(eostype);
        const double epsilon = internal::epsilon(eostype);
        const double sigma = internal::sigma(eostype);
        const auto alpha = internal::alpha(eostype);

        // Calculate the parameters `a` of the cubic equation of state for each species
        ThermoVector a(nspecies);
        ThermoVector aT(nspecies);
        ThermoVector aTT(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tc = critical_temperatures[i];
            const double Pc = critical_pressures[i];
            const double omega = acentric_factors[i];
            const double factor = Psi*R*R*(Tc*Tc)/(Pc*Pc);
            const ThermoScalar Tr = T/Tc;
            const ThermoScalar alpha_val, alpha_ddt, alpha_d2dt2;
            std::tie(alpha_val, alpha_ddt, alpha_d2dt2) = alpha(Tr, omega);
            a[i] = factor * alpha_val;
            aT[i] = factor * alpha_ddt;
            aTT[i] = factor * alpha_d2dt2;
        };

        // Calculate the parameters `b` of the cubic equation of state for each species
        ThermoVector b(nspecies);
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
        ChemicalScalar amix(nspecies);
        ChemicalScalar amixT(nspecies);
        ChemicalScalar amixTT(nspecies);
        ChemicalVector abar(nspecies, nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            for(unsigned j = 0; j < nspecies; ++j)
            {
                const ThermoScalar r = kres.k.empty() ? 1.0 : 1.0 - kres.k[i][j];
                const ThermoScalar rT = kres.kT.empty() ? 0.0 : -kres.kT[i][j];
                const ThermoScalar rTT = kres.kTT.empty() ? 0.0 : -kres.kTT[i][j];

                const ThermoScalar s = sqrt(a[i]*a[j]);
                const ThermoScalar sT = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
                const ThermoScalar sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;

                const ThermoScalar aij = r*s;
                const ThermoScalar aijT = rT*s + r*sT;
                const ThermoScalar aijTT = rTT*s + 2*rT*sT + r*sTT;

                amix += x[i] * x[j] * aij;
                amixT += x[i] * x[j] * aijT;
                amixTT += x[i] * x[j] * aijTT;

                abar[i] += 2 * x[j] * aij;
            }
        }

        // Finalize the calculation of `abar`
        for(unsigned i = 0; i < nspecies; ++i)
            abar[i] -= amix;

        // Calculate the parameter `bmix` of the cubic equation of state
        ChemicalScalar bmix(nspecies);
        ThermoVector bbar(nspecies);
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
        const ChemicalScalar qTT = qT*qT/q + q*(1/(T*T) + aTT/a - aT*aT/(a*a));

        // Calculate the coefficients A, B, C of the cubic equation of state
        const ChemicalScalar A = (epsilon + sigma - 1)*beta - 1;
        const ChemicalScalar B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon - sigma + q)*beta;
        const ChemicalScalar C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

        // Calculate the partial temperature derivative of the coefficients A, B, C
        const ChemicalScalar AT = (epsilon + sigma - 1)*betaT;
        const ChemicalScalar BT = 2*(epsilon*sigma - epsilon - sigma)*beta*betaT - qT*beta - (epsilon - sigma + q)*betaT;
        const ChemicalScalar CT = -3*epsilon*sigma*beta*beta*betaT - qT*beta*beta - 2*(epsilon*sigma + q)*beta*betaT;

        // Determine the appropriate initial guess for the cubic equation of state
        const double Z0 = isvapor ? 1.0 : beta.val;

        // Define the non-linear function and its derivative
        const auto f = [&](double Z) -> double { return Z*Z*Z + A.val*Z*Z + B.val*Z + C.val; };
        const auto g = [&](double Z) -> double { return 3*Z*Z + 2*A.val*Z + B.val; };
        const auto tol = 1e-6;
        const auto maxiter = 10;

        // Calculate the compressibility factor Z
        ChemicalScalar Z(nspecies);
        Z.val = newton(f, g, Z0, tol, maxiter);

        // Calculate the partial derivatives of Z (dZdT, dZdP, dZdn)
        const double factor = -1.0/(3*Z.val*Z.val + 2*A.val*Z.val + B.val);
        Z.ddt = factor * (A.ddt*Z.val*Z.val + B.ddt*Z.val + C.ddt);
        Z.ddp = factor * (A.ddp*Z.val*Z.val + B.ddp*Z.val + C.ddp);
        for(unsigned i = 0; i < nspecies; ++i)
            Z.ddn[i] = factor * (A.ddn[i]*Z.val*Z.val + B.ddn[i]*Z.val + C.ddn[i]);

        // Calculate the partial temperature derivative of Z
        const ChemicalScalar ZT = -(AT*Z*Z + BT*Z + CT)/(3*Z*Z + 2*A*Z + B);

        // Calculate the integration factor I and its temperature derivative IT
        const ChemicalScalar I = (epsilon != sigma) ?
            log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon) :
                beta/(Z + epsilon*beta);

        const ChemicalScalar IT = (epsilon != sigma) ?
            ((ZT + sigma*betaT)/(Z + sigma*beta) - (ZT + epsilon*betaT)/(Z + epsilon*beta))/(sigma - epsilon) :
                I*(betaT/beta - (ZT + epsilon*betaT)/(Z + epsilon*beta));

        Result result;
        ChemicalScalar& V = result.molar_volume;
        ChemicalScalar& G_res = result.residual_molar_gibbs_energy;
        ChemicalScalar& H_res = result.residual_molar_enthalpy;
        ChemicalScalar& Cp_res = result.residual_molar_heat_capacity_cp;
        ChemicalScalar& Cv_res = result.residual_molar_heat_capacity_cv;

        const ChemicalScalar dPdT = P*(1/T + ZT/Z);
        const ChemicalScalar dVdT = V*(1/T + ZT/Z);

        // Calculate the partial molar Zi for each species
        V = Z*R*T/P;
        G_res = R*T*(Z - 1 - log(Z - beta) - q*I);
        H_res = R*T*(Z - 1 + T*qT*I);
        Cp_res = R*T*(ZT + qT*I + T*qTT + T*qT*IT) + H_res/T;
        Cv_res = Cp_res - T*dPdT*dVdT + R;

        for(unsigned i = 0; i < nspecies; ++i)
        {
            const ChemicalScalar ai = abar[i];
            const ChemicalScalar bi = bbar[i];
            const ChemicalScalar betai = P*bi/(R*T);
            const ChemicalScalar qi = q*(1 + ai/amix + bi/bmix);
            const ChemicalScalar Ai = (epsilon + sigma - 1)*betai - 1;
            const ChemicalScalar Bi = (epsilon*sigma - epsilon - sigma)*betai*betai - (epsilon - sigma + qi)*betai;
            const ChemicalScalar Ci = -epsilon*sigma*betai*betai*betai - (epsilon*sigma + qi)*betai*betai;
            const ChemicalScalar Zi = -(Ai*Z*Z + (Bi + B)*Z + Ci + 2*C)/(3*Z*Z + 2*A*Z + B);
            const ChemicalScalar Ii = (epsilon != sigma) ?
                I + ((Zi + sigma*betai)/(Z + sigma*beta) - (Zi + epsilon*betai)/(Z + epsilon*beta))/(sigma - epsilon) :
                I * (1 + betai/beta - (Zi + epsilon*betai)/(Z + epsilon*beta));

            result.fugacity_coefficients[i] = Zi - (Zi - betai)/(Z - beta) - log(Z - beta) - qi*I - q*Ii + q*I;

        }

    }
};

CubicEOS::CubicEOS(unsigned nspecies)
: pimpl(new Impl(nspecies))
{}

auto CubicEOS::numSpecies() const -> unsigned
{
    return pimpl->nspecies;
}

auto CubicEOS::setPhaseIsLiquid() -> void
{
    pimpl->isvapor = false;
}

auto CubicEOS::setPhaseIsVapor() -> void
{
    pimpl->isvapor = true;
}

auto CubicEOS::setCriticalTemperatures(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the critical "
        "temperatures of the species in CubicEOS.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " were given.");

    pimpl->critical_temperatures = values;
}

auto CubicEOS::setCriticalPressures(const std::vector<double>& values) -> void
{
    Assert(values.size() == numSpecies(), "Cannot set the critical "
        "pressures of the species in CubicEOS.", "Expecting " +
        std::to_string(numSpecies()) + " values, but only " +
        std::to_string(values.size()) + " were given.");

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

auto CubicEOS::setType(Type type) -> void
{
    pimpl->eostype = type;
}

//auto CubicEOS::addBinaryInteractionParams(Index i, Index j, const std::vector<double>& kij) -> void
//{
//    pimpl->ij.push_back(std::make_tuple(i, j));
//    pimpl->kij.push_back(kij);
//}

auto CubicEOS::operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> Result
{
    return pimpl->operator()(T, P, x);
}


// todo remove this
//ChemicalVector n = ChemicalVector::Composition(nvec);
//ChemicalScalar nt = sum(n);
//
//x = n/nt;
//
//nt.val = sum(n.val);
//nt.ddt = sum(n.ddt);
//nt.ddp = sum(n.ddp);
//nt.ddn = sumcols(n.ddn);

} // namespace Reaktoro
