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

    /// Calculate the thermodynamic properties of the phase.
    /// @param T The temperature of the phase (in units of K)
    /// @param P The pressure of the phase (in units of Pa)
    /// @param x The molar fractions of the species of the phase (in units of mol/mol)
    auto operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> EOSResult
    {
        // The table of binary interaction parameters
        Table2D<ThermoScalar> k = calculate_k(T);

        // The parameters `a` of the cubic equation of state for each species
        ThermoVector a = calculate_a(T);

        // The parameter `amix` of the cubic equation of state
        ChemicalScalar amix = calculate_amix(x, a, k);

        // The parameter `bmix` of the cubic equation of state
        ChemicalScalar bmix = calculate_bmix(x);


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
