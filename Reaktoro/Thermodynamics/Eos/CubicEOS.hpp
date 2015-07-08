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
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

/// Defines a cubic equation of state and calculates thermodynamic properties of a fluid phase.
class CubicEOS
{
public:
    /// Defines the types of supported cubic equation of states.
    enum Type
    {
        VanDerWaals, RedlichKwong, SoaveRedlichKwong, PengRobinson,
    };

    struct InteractionParamsResult
    {
        Table2D<ThermoScalar> k;

        Table2D<ThermoScalar> kT;

        Table2D<ThermoScalar> kTT;
    };

    struct InteractionParamsArgs
    {
        const ChemicalScalar& T;

        const ThermoVector& a;

        const ThermoVector& aT;

        const ThermoVector& aTT;

        const ThermoVector& b;
    };

    using InteractionParamsFunction =
        std::function<InteractionParamsResult(const InteractionParamsArgs&)>;

    struct Result
    {
        /// The molar volume of the phase (in units of m3/mol).
        ChemicalScalar molar_volume;

        /// The residual molar Gibbs energy of the phase (in units of J/mol).
        ChemicalScalar residual_molar_gibbs_energy;

        /// The residual molar enthalpy of the phase (in units of J/mol).
        ChemicalScalar residual_molar_enthalpy;

        /// The residual molar heat capacity at constant pressure of the phase (in units of J/(mol*K)).
        ChemicalScalar residual_molar_heat_capacity_cp;

        /// The residual molar heat capacity at constant volume of the phase (in units of J/(mol*K)).
        ChemicalScalar residual_molar_heat_capacity_cv;

        /// The fugacity coefficients of the species in the phase.
        ChemicalVector fugacity_coefficients;
    };

    /// Construct a CubicEOS instance with given number of species.
    /// @param nspecies The number of species in the phase.
    explicit CubicEOS(unsigned nspecies);

    /// Return the number of species in the phase.
    auto numSpecies() const -> unsigned;

    /// Set the equation of state to compute properties for a liquid phase.
    auto setPhaseIsLiquid() -> void;

    /// Set the equation of state to compute properties for a vapor phase.
    auto setPhaseIsVapor() -> void;

    /// Set the critical temperatures of the species (in units of K).
    auto setCriticalTemperatures(const std::vector<double>& values) -> void;

    /// Set the critical pressures of the species (in units of Pa).
    auto setCriticalPressures(const std::vector<double>& values) -> void;

    /// Set the acentric factors of the species.
    auto setAcentricFactors(const std::vector<double>& values) -> void;

    /// Set the function that calculates the interaction parameters kij and its temperature derivatives.
    /// @see InteractionParamFunction, InteractionParamArgs, InteractionParamResult
    auto setInteractionParamsFunction(const InteractionParamsFunction& func) -> void;

    /// Set the type of the cubic equation of state (default: PengRobinson).
    /// @see Type
    auto setType(Type type) -> void;

    /// Add binary interaction parameters for the the species pair *(i, j)*.
    /// @param i The index of the species *i*
    /// @param j The index of the species *j*
    /// @param kij The coefficients of the binary interaction parameter equation \f$ k_{ij}(T) \f$
    ///        for the pair of species *i* and *j*
    auto addBinaryInteractionParams(Index i, Index j, const std::vector<double>& kij);

    /// Calculate the thermodynamic properties of the phase.
    /// @param T The temperature of the phase (in units of K)
    /// @param P The pressure of the phase (in units of Pa)
    /// @param x The molar fractions of the species of the phase (in units of mol/mol)
    auto operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> EOSResult;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
