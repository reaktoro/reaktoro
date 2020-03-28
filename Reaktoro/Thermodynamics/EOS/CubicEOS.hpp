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

#pragma once

// C++ includes
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Common/TableUtils.hpp>
#include <Reaktoro/Thermodynamics/EOS/PhaseIdentification.hpp>

namespace Reaktoro {

/// Defines a cubic equation of state and calculates thermodynamic properties of a fluid phase.
class CubicEOS
{
public:
    /// Defines the enumeration of available cubic EOS models.
    enum Model
    {
        VanDerWaals, RedlichKwong, SoaveRedlichKwong, PengRobinson,
    };

    /// Class to define or store Binary Interaction Parameters (BIPs) from a calculation or input.
    /// Note that the BIPs can depend on the temperature, thus kT and kTT should be also provided.
    struct InteractionParamsResult
    {
        MatrixXd k;

        MatrixXd kT;

        MatrixXd kTT;
    };

    /// Function wrapper to calculate (temperature-dependent) binary interaction parameters.
    using InteractionParamsFunction =
        std::function<InteractionParamsResult(const double&)>;

    /// Parameters to be passed to the Cubic Equation of State
    struct Params
    {
        Model model = PengRobinson;

        /// If both Gaseous and Liquid phases are in the system, it is recommended to configure a
        /// robust phase identification method such as GibbsEnergyAndEquationOfStateMethod for BOTH
        /// phases.
        PhaseIdentificationMethod phase_identification_method = PhaseIdentificationMethod::None;
    };

    struct InteractionParamsResult
    {
        Table2D<real> k;

        Table2D<real> kT;

        Table2D<real> kTT;
    };

    struct InteractionParamsArgs
    {
        const real& T;

        const VectorXr& a;

        const VectorXr& aT;

        const VectorXr& aTT;

        VectorXrConstRef b;
    };

    using InteractionParamsFunction =
        std::function<InteractionParamsResult(const InteractionParamsArgs&)>;

    struct Result
    {
        /// Construct a default Result instance
        Result();

        /// Construct a Result instance with zero entries
        /// @param nspecies The number of species
        explicit Result(unsigned nspecies);

        /// The molar volume of the phase (in units of m3/mol).
        real molar_volume;

        /// The residual molar Gibbs energy of the phase (in units of J/mol).
        real residual_molar_gibbs_energy;

        /// The residual molar enthalpy of the phase (in units of J/mol).
        real residual_molar_enthalpy;

        /// The residual molar heat capacity at constant pressure of the phase (in units of J/(mol*K)).
        real residual_molar_heat_capacity_cp;

        /// The residual molar heat capacity at constant volume of the phase (in units of J/(mol*K)).
        real residual_molar_heat_capacity_cv;

        /// The partial molar volumes of the species in the phase (in units of m3/mol).
        VectorXd partial_molar_volumes;

        /// The residual partial molar Gibbs energies of the species in the phase (in units of J/mol).
        VectorXd residual_partial_molar_gibbs_energies;

        /// The residual partial molar enthalpies of the species in the phase (in units of J/mol).
        VectorXd residual_partial_molar_enthalpies;

        /// The fugacity coefficients of the species in the phase.
        VectorXd ln_fugacity_coefficients;
    };

    /// Construct a CubicEOS instance with given number of species.
    /// @param nspecies The number of species in the phase.
    explicit CubicEOS(unsigned nspecies, Params params);

    /// Construct a copy of a CubicEOS instance
    CubicEOS(const CubicEOS& other);

    /// Destroy this CubicEOS instance
    virtual ~CubicEOS();

    /// Assign a CubicEOS instance to this
    auto operator=(CubicEOS other) -> CubicEOS&;

    /// Return the number of species in the phase.
    auto numSpecies() const -> unsigned;

    /// Set the type of the cubic equation of state (default: PengRobinson).
    /// @see Model
    auto setModel(Model model) -> void;

    /// Set the equation of state to compute properties for a liquid phase.
    auto setPhaseAsLiquid() -> void;

    /// Set the equation of state to compute properties for a vapor phase.
    auto setPhaseAsVapor() -> void;

    /// Set the critical temperatures of the species (in units of K).
    auto setCriticalTemperatures(const std::vector<double>& values) -> void;

    /// Set the critical pressures of the species (in units of Pa).
    auto setCriticalPressures(const std::vector<double>& values) -> void;

    /// Set the acentric factors of the species.
    auto setAcentricFactors(const std::vector<double>& values) -> void;

    /// Set the function that calculates the interaction parameters kij (k) 
    /// and its temperature derivatives (kT and kTT).
    /// @param func A wrapper function to calculate the BIPs matrices
    /// @see InteractionParamFunction, InteractionParamArgs, InteractionParamResult
    auto setInteractionParamsFunction(const InteractionParamsFunction& func) -> void;

    /// Calculate the thermodynamic properties of the phase.
    /// @param T The temperature of the phase (in units of K)
    /// @param P The pressure of the phase (in units of Pa)
    /// @param x The mole fractions of the species in the phase (in units of mol/mol)
    auto operator()(const real& T, const real& P, const VectorXd& x) -> Result;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
