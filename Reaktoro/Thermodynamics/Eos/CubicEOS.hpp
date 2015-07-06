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

#pragma once

// C++ includes
#include <string>
#include <vector>
#include <tuple>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

struct EosResult
{
    /// The molar volume of the phase (in units of m3/mol).
    ChemicalScalar molar_volume;

    /// The fugacity coefficients of the species in the phase.
    ChemicalVector fugacity_coefficients;

    /// The residual molar Gibbs energy of the phase (in units of J/mol).
    ChemicalScalar residual_molar_gibbs_energy;

    /// The residual molar enthalpy energy of the phase (in units of J/mol).
    ChemicalScalar residual_molar_enthalpy_energy;
};

enum class CubicEosType
{
    VanDerWaals,
    RedlichKwong,
    SoaveRedlichKwong,
    PengRobinson,
};

struct CubicEosParams
{
    /// The type that defines a table of binary interaction parameters for pair of species.
    using TableIndexIndexParams =
        std::vector<std::tuple<Index, Index, std::vector<double>>>;

    /// The type of the phase, which is either `'v'` for vapor, or `'l'` for liquid.
    char phasetype = 'v';

    /// The type of the cubic equation of state.
    CubicEosType eostype = CubicEosType::PengRobinson;

    /// The critical temperatures of the species (in units of K).
    std::vector<double> critical_temperatures;

    /// The critical pressures of the species (in units of Pa).
    std::vector<double> critical_pressures;

    /// The acentric factors of the species.
    std::vector<double> acentric_factors;

    /// The table of binary interaction parameters.
    /// The i-th entry in this table gives a tuple `(index1, index2, params)`,
    /// where `index1` and `index2` are the indices of the species for which `params`
    /// hold their binary interaction parameters. The number of interaction parameters
    /// per pair of species are either **one** or **three**. If three parameters are given,
    /// these are used in the following correlation equation for temperature:
    /// TODO add equation here.
    TableIndexIndexParams binary_interaction_params;
};

/// Defines a cubic equation of state and calculates thermodynamic properties of a fluid phase.
class CubicEos
{
    /// Construct a default CubicEos instance.
    CubicEos();

    /// Construct a custom CubicEos instance with given parameters.
    CubicEos(const CubicEosParams& params);

    /// Calculate the thermodynamic properties of the phase.
    /// @param T The temperature of the phase (in units of K)
    /// @param P The pressure of the phase (in units of Pa)
    /// @param x The molar fractions of the species of the phase (in units of mol/mol)
    auto operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> EosResult;
};

namespace internal {

auto alpha(CubicEosType type) -> std::function<ThermoScalar(const ThermoScalar&, double)>
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
        case CubicEosType::VanDerWaals: return alphaVDW;
        case CubicEosType::RedlichKwong: return alphaRK;
        case CubicEosType::SoaveRedlichKwong: return alphaSRK;
        case CubicEosType::PengRobinson: return alphaPR;
        default: return alphaPR;
    }
}

auto sigma(CubicEosType type) -> double
{
    switch(type)
    {
        case CubicEosType::VanDerWaals: return 0.0;
        case CubicEosType::RedlichKwong: return 1.0;
        case CubicEosType::SoaveRedlichKwong: return 1.0;
        case CubicEosType::PengRobinson: return 1.0 + 1.41421356237;
        default: return 1.0 + 1.41421356237;
    }
}

auto epsilon(CubicEosType type) -> double
{
    switch(type)
    {
        case CubicEosType::VanDerWaals: return 0.0;
        case CubicEosType::RedlichKwong: return 0.0;
        case CubicEosType::SoaveRedlichKwong: return 0.0;
        case CubicEosType::PengRobinson: return 1.0 - 1.41421356237;
        default: return 1.0 - 1.41421356237;
    }
}

auto Omega(CubicEosType type) -> double
{
    switch(type)
    {
        case CubicEosType::VanDerWaals: return 1.0/8.0;
        case CubicEosType::RedlichKwong: return 0.08664;
        case CubicEosType::SoaveRedlichKwong: return 0.08664;
        case CubicEosType::PengRobinson: return 0.07780;
        default: return 0.07780;
    }
}

auto Psi(CubicEosType type) -> double
{
    switch(type)
    {
        case CubicEosType::VanDerWaals: return 27.0/64.0;
        case CubicEosType::RedlichKwong: return 0.42748;
        case CubicEosType::SoaveRedlichKwong: return 0.42748;
        case CubicEosType::PengRobinson: return 0.45724;
        default: return 0.45724;
    }
}

} // namespace internal

auto a = [](const ThermoScalar& T, const CubicEosParams& params) -> ThermoVector
{
    const unsigned nspecies = params.nspecies;
    const double R = universalGasConstant;
    const double R2 = R*R;
    const double Psi = internal::Psi(params.eostype);
    const auto alpha = internal::alpha(params.eostype);
    ThermoVector res(nspecies);
    ThermoScalar Tr;
    for(unsigned i = 0; i < nspecies; ++i)
    {
        const double Tc = params.critical_temperatures[i];
        const double Pc = params.critical_pressures[i];
        const double omega = params.acentric_factors[i];
        Tr = T/Tc;
        res[i] = Psi * alpha(Tr, omega) * R2 * (Tc * Tc)/(Pc * Pc);
    };
    return res;
};

auto amix = [](const ChemicalVector& x, const ThermoVector& a, const std::vector<std::vector<ThermoScalar>>& k) -> ChemicalScalar
{
    const unsigned nspecies = size(x);
    ChemicalScalar res(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
    {
        for(unsigned j = 0; j < nspecies; ++j)
        {
            ThermoScalar kij = k[i][j];
            ThermoScalar aij = (1 - kij) * sqrt(a[i] * a[j]);
            res += x[i] * x[j] * aij;
        }
    }
    return res;
};

auto bmix = [](const ChemicalVector& x, const Vector& b) -> ChemicalScalar
{
    const unsigned nspecies = size(x);
    ChemicalScalar res(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        res += x[i] * b[i];
    return res;
};

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
