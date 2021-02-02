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
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/TableUtils.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Math/Matrix.hpp>
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
        /// The BIPs matrix. The size must be (n, n), where n is the number of species
        MatrixXd k;

        /// The derivative of each k entry w.r.t T.
        MatrixXd kT;

        /// The derivative of each kT entry w.r.t T.
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

        /// Wrapper function to calculate binary interaction parameters for mixture rules using a 
        /// Cubic EOS. This function should have the temperature T as input 
        /// and should return a InteractionParamsResult.
        /// @see InteractionParamsResult
        InteractionParamsFunction binary_interaction_values;
    };

    struct Result
    {
        /// Construct a default Result instance
        Result();

        /// Construct a Result instance with zero entries
        /// @param nspecies The number of species
        explicit Result(unsigned nspecies);

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

        /// The partial molar volumes of the species in the phase (in units of m3/mol).
        ChemicalVector partial_molar_volumes;

        /// The residual partial molar Gibbs energies of the species in the phase (in units of J/mol).
        ChemicalVector residual_partial_molar_gibbs_energies;

        /// The residual partial molar enthalpies of the species in the phase (in units of J/mol).
        ChemicalVector residual_partial_molar_enthalpies;

        /// The fugacity coefficients of the species in the phase.
        ChemicalVector ln_fugacity_coefficients;
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

    /// Sanity check to verify if BIPs matrices have proper dimensions. Considering that the phase has
    /// n species, the BIP matricies k, kT and kTT should have (n, n) as dimensions.
    /// @see setInteractionParamsFunction
    auto sanityCheckInteractionParamsFunction(const InteractionParamsFunction& func) -> void;

    /// Calculate the thermodynamic properties of the phase.
    /// @param T The temperature of the phase (in units of K)
    /// @param P The pressure of the phase (in units of Pa)
    /// @param x The mole fractions of the species in the phase (in units of mol/mol)
    auto operator()(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& x) -> Result;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
