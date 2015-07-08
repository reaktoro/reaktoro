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

auto alpha(CubicEOS::Type type) -> std::function<ThermoScalar(const ThermoScalar&, double)>
{
    auto alphaVDW = [](const ThermoScalar& Tr, double omega) -> ThermoScalar
    {
        return 1.0;
    };

    auto alphaRK = [](const ThermoScalar& Tr, double omega) -> ThermoScalar
    {
        return 1.0/sqrt(Tr);
    };

    auto alphaSRK = [](const ThermoScalar& Tr, double omega) -> ThermoScalar
    {
        double factor = 0.480 + 1.574*omega - 0.176*omega*omega;
        ThermoScalar aux = 1.0 + factor*(1.0 - sqrt(Tr));
        return aux*aux;
    };

    auto alphaPR = [](const ThermoScalar& Tr, double omega) -> ThermoScalar
    {
        double factor = 0.37464 + 1.5422*omega - 0.26992*omega*omega;
        ThermoScalar aux = 1.0 + factor*(1.0 - sqrt(Tr));
        return aux*aux;
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
        case CubicEOS::PengRobinson: return 0.07780;
        default: return 0.07780;
    }
}

auto Psi(CubicEOS::Type type) -> double
{
    switch(type)
    {
        case CubicEOS::VanDerWaals: return 27.0/64.0;
        case CubicEOS::RedlichKwong: return 0.42748;
        case CubicEOS::SoaveRedlichKwong: return 0.42748;
        case CubicEOS::PengRobinson: return 0.45724;
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

    /// The pairs of species indices for the binary interaction parameters `kij`
    std::vector<std::tuple<Index, Index>> ij;

    /// The binary intectaction parameters `kij`. The number of interaction parameters
    /// per pair of species are either **one** or **three**. If three parameters are given,
    /// these are used in the following correlation equation for temperature:
    std::vector<std::vector<double>> kij;

    /// Construct a CubicEOS::Impl instance.
    Impl(unsigned nspecies)
    : nspecies(nspecies)
    {}

    /// Calculate the table of binary interaction parameters.
    auto calculate_k(const ThermoScalar& T) -> Table2D<ThermoScalar>
    {
        Table2D<ThermoScalar> k = table2D<ThermoScalar>(nspecies, nspecies);
        for(unsigned ipair = 0; ipair < ij.size(); ++ipair)
        {
            const Index i = std::get<0>(ij[ipair]);
            const Index j = std::get<1>(ij[ipair]);
            const double kval = kij[ipair][0];
            k[i][j].val = kval;
        }
        return k;
    };

    /// Calculate the parameter `a` for each species in the phase.
    auto calculate_a(const ThermoScalar& T) -> ThermoVector
    {
        ThermoVector a(nspecies);
        const double R = universalGasConstant;
        const double R2 = R*R;
        const double Psi = internal::Psi(eostype);
        const auto alpha = internal::alpha(eostype);
        ThermoScalar Tr;
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tc = critical_temperatures[i];
            const double Pc = critical_pressures[i];
            const double omega = acentric_factors[i];
            Tr = T/Tc;
            a[i] = Psi * alpha(Tr, omega) * R2 * (Tc * Tc)/(Pc * Pc);
        };
        return a;
    };

    // Calculate the parameter `amix` of the cubic equation of state.
    auto calculate_amix(const ChemicalVector& x, const ThermoVector& a, const Table2D<ThermoScalar>& k) -> ChemicalScalar
    {
        ChemicalScalar amix(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            for(unsigned j = 0; j < nspecies; ++j)
            {
                ThermoScalar kij = k[i][j];
                ThermoScalar aij = (1 - kij) * sqrt(a[i] * a[j]);
                amix += x[i] * x[j] * aij;
            }
        }
        return amix;
    };

    /// Calculate the parameter `b` for each species in the phase.
    auto calculate_bmix(const ChemicalVector& x) -> ChemicalScalar
    {
        ChemicalScalar bmix(nspecies);
        const double R = universalGasConstant;
        const double Omega = internal::Omega(eostype);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tci = critical_temperatures[i];
            const double Pci = critical_pressures[i];
            const double bi = Omega*R*Tci/Pci;
            bmix += x[i] * bi;
        }
        return bmix;
    };

    auto calculate_Z(const ChemicalScalar& beta, const ChemicalScalar& q) -> ChemicalScalar
    {
        const double epsilon = internal::epsilon(eostype);
        const double sigma = internal::sigma(eostype);
        const ChemicalScalar A = (epsilon + sigma - 1)*beta - 1;
        const ChemicalScalar B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon - sigma + q)*beta;
        const ChemicalScalar C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

        const double Z0 = isvapor ? 1.0 : beta.val;

        // Calculate the
        const auto f = [](double Z) -> double { return Z*Z*Z + A.val*Z*Z + B.val*Z + C.val; };
        const auto g = [](double Z) -> double { return 3*Z*Z + 2*A.val*Z + B.val; };
        const auto tol = 1e-6;
        const auto maxiter = 10;
        ChemicalScalar Z(nspecies);
        Z.val = newton(f, g, Z0, tol, maxiter);

        const double Z2 = Z.val*Z.val;
        const double factor = -1.0/(3*Z2 + 2*A.val*Z.val + B.val);

        Z.ddt = factor * (A.ddt*Z2 + B.ddt*Z.val + C.ddt);
        Z.ddp = factor * (A.ddp*Z2 + B.ddp*Z.val + C.ddp);
        for(unsigned i = 0; i < nspecies; ++i)
            Z.ddn[i] = factor * (A.ddn[i]*Z2 + B.ddn[i]*Z.val + C.ddn[i]);

        return Z;
    };

    auto calculate_V(const ThermoScalar& T, const ThermoScalar& P, const ChemicalScalar& Z) -> ChemicalScalar
    {
        const double R = universalGasConstant;
        return R*T*Z/P;
    }

    auto operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> EOSResult
    {
        // Auxiliary variables
        const double R = universalGasConstant;
        const double Psi = internal::Psi(eostype);
        const double Omega = internal::Omega(eostype);
        const double epsilon = internal::epsilon(eostype);
        const double sigma = internal::sigma(eostype);
        const auto alpha = internal::alpha(eostype);

        // Initialize the table of binary interaction parameters
        Table2D<ThermoScalar> k = table2D<ThermoScalar>(nspecies, nspecies);
        for(unsigned ipair = 0; ipair < ij.size(); ++ipair)
        {
            const Index i = std::get<0>(ij[ipair]);
            const Index j = std::get<1>(ij[ipair]);
            const double kval = kij[ipair][0];
            k[i][j].val = kval;
        }

        // Calculate the parameters `a` of the cubic equation of state for each species
        ThermoVector a(nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            const double Tc = critical_temperatures[i];
            const double Pc = critical_pressures[i];
            const double omega = acentric_factors[i];
            const ThermoScalar Tr = T/Tc;
            a[i] = Psi * alpha(Tr, omega) * R*R * (Tc * Tc)/(Pc * Pc);
        };

        // Calculate the parameter `amix` of the phase and the partial molar parameters `abar` of each species
        ChemicalScalar amix(nspecies);
        ChemicalVector abar(nspecies, nspecies);
        for(unsigned i = 0; i < nspecies; ++i)
        {
            for(unsigned j = 0; j < nspecies; ++j)
            {
                const ThermoScalar kij = k[i][j];
                const ThermoScalar aij = (1 - kij) * sqrt(a[i] * a[j]);
                amix += x[i] * x[j] * aij;
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

        // Calculate auxiliary quantities `beta` and `q`
        const ChemicalScalar beta = P*bmix/(R*T);
        const ChemicalScalar q = amix/(bmix*R*T);

        // Calculate the coefficients of the cubic equation of state
        const ChemicalScalar A = (epsilon + sigma - 1)*beta - 1;
        const ChemicalScalar B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon - sigma + q)*beta;
        const ChemicalScalar C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

        // Determine the appropriate initial guess for the cubic equation of state
        const double Z0 = isvapor ? 1.0 : beta.val;

        // Define the non-linear function and its derivative
        const auto f = [](double Z) -> double { return Z*Z*Z + A.val*Z*Z + B.val*Z + C.val; };
        const auto g = [](double Z) -> double { return 3*Z*Z + 2*A.val*Z + B.val; };
        const auto tol = 1e-6;
        const auto maxiter = 10;

        // Calculate the compressibility factor Z
        ChemicalScalar Z(nspecies);
        Z.val = newton(f, g, Z0, tol, maxiter);

        // Calculate the derivatives of Z w.r.t. (T, P, n)
        const double factor = -1.0/(3*Z.val*Z.val + 2*A.val*Z.val + B.val);
        Z.ddt = factor * (A.ddt*Z.val*Z.val + B.ddt*Z.val + C.ddt);
        Z.ddp = factor * (A.ddp*Z.val*Z.val + B.ddp*Z.val + C.ddp);
        for(unsigned i = 0; i < nspecies; ++i)
            Z.ddn[i] = factor * (A.ddn[i]*Z.val*Z.val + B.ddn[i]*Z.val + C.ddn[i]);

        // Calculate the integration factor I
        ChemicalScalar I = (epsilon != sigma) ?
            log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon) :
                beta/(Z + epsilon*beta);

        // Calculate the partial molar Zi for each species
        Result result;
        result.molar_volume = Z*R*T/P;
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
                I * (1 + betai/beta + (Zi + epsilon*betai)/(Z + epsilon*beta));

            result.residual_partial_molar_volumes[i] = Zi*R*T/P;
            result.residual_partial_molar_volumes[i] = Zi*R*T/P;
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

auto CubicEOS::setType(Type type) -> void
{
    pimpl->eostype = type;
}

auto CubicEOS::addBinaryInteractionParams(Index i, Index j, const std::vector<double>& kij)
{
    pimpl->ij.push_back(std::make_tuple(i, j));
    pimpl->kij_params.push_back(kij);
}

auto CubicEOS::operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> EOSResult
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
