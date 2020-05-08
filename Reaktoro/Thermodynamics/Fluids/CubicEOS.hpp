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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Thermodynamics/Fluids/PhaseIdentification.hpp>

namespace Reaktoro {

/// The properties calculated by evaluating a cubic equation of state.
struct CubicEOSProps
{
    /// The residual molar volume of the phase (in m3/mol).
    real Vres = {};

    /// The temperature derivative of the residual molar volume at constant pressure (in m3/(mol*K)).
    real VresT = {};

    /// The pressure derivative of the residual molar volume constant temperature (in m3/(mol*Pa)).
    real VresP = {};

    /// The residual molar Gibbs energy of the phase (in J/mol).
    real Gres = {};

    /// The residual molar enthalpy of the phase (in J/mol).
    real Hres = {};

    /// The residual molar heat capacity at constant pressure of the phase (in J/(mol*K)).
    real Cpres = {};

    /// The residual molar heat capacity at constant volume of the phase (in J/(mol*K)).
    real Cvres = {};

    /// The ln fugacity coefficients of the species in the phase.
    ArrayXr ln_phi;
};

/// The options for the cubic equation of state models.
enum class CubicEOSModel
{
    VanDerWaals,
    RedlichKwong,
    SoaveRedlichKwong,
    PengRobinson
};

/// The interaction parameters for the cubic equation of state.
struct CubicEOSInteractionParams
{
    MatrixXr k;

    MatrixXr kT;

    MatrixXr kTT;
};

/// The arguments needed in the function that evaluates interaction parameters for the cubic equation of state.
struct CubicEOSInteractionParamsArgs
{
    const real& T;

    ArrayXrConstRef a;

    ArrayXrConstRef aT;

    ArrayXrConstRef aTT;

    ArrayXrConstRef b;
};

/// The function that evaluates interaction parameters for the cubic equation of state model.
using CubicEOSInteractionParamsFn = std::function<CubicEOSInteractionParams(const CubicEOSInteractionParamsArgs&)>;

/// Calculates thermodynamic properties of fluid phases based on a cubic equation of state model.
class CubicEOS
{
public:
    /// The arguments needed to construct a CubicEOS object.
    struct Args
    {
        /// The number of species in the fluid phase.
        const Index nspecies = {};

        /// The critical temperatures of the substances (in K).
        ArrayXrConstRef Tcr;

        /// The critical pressures of the substances (in Pa).
        ArrayXrConstRef Pcr;

        /// The acentric factor of the substances.
        ArrayXrConstRef omega;

        /// The fluid type for which the equation of state should be confifured.
        CubicEOSFluidType fluidtype = CubicEOSFluidType::Vapor;

        /// The cubic equation of state model to be used.
        CubicEOSModel model = CubicEOSModel::PengRobinson;

        /// The function that calculates interaction parameters @eq{k_{ij}} in @eq{a_{ij}=(1-k_{ij})(a_{i}a_{j})^{1/2}}.
        CubicEOSInteractionParamsFn interaction_params_fn;

        /// The method to identify whether liquid or vapor phases is stable.
        /// If both vapor and liquid phases are considered in the system, it is
        /// recommended to configure a robust phase identification method such as
        /// GibbsEnergyAndEquationOfStateMethod for BOTH phases.
        PhaseIdentificationMethod phase_identification_method = PhaseIdentificationMethod::None;
    };

    /// Construct a CubicEOS instance with given number of species.
    /// @param nspecies The number of species in the phase.
    explicit CubicEOS(const Args& args);

    /// Construct a copy of a CubicEOS instance
    CubicEOS(const CubicEOS& other);

    /// Destroy this CubicEOS instance
    virtual ~CubicEOS();

    /// Assign a CubicEOS instance to this
    auto operator=(CubicEOS other) -> CubicEOS&;

    /// Return the number of species in the phase.
    auto numSpecies() const -> unsigned;

    /// Set the cubic equation of state model.
    auto setModel(CubicEOSModel model) -> void;

    /// Set the fluid type, liquid or vapor, for which the equation of state should be confifured.
    auto setFluidType(CubicEOSFluidType fluidtype) -> void;

    /// Set the function that calculates interaction parameters @eq{k_{ij}} in @eq{a_{ij}=(1-k_{ij})(a_{i}a_{j})^{1/2}}.
    /// @see CubicEOSInteractionParamsFn, CubicEOSInteractionParams
    auto setInteractionParamsFunction(const CubicEOSInteractionParamsFn& func) -> void;

    /// Set the method to identify whether liquid or vapor phases is stable.
    auto setStablePhaseIdentificationMethod(const PhaseIdentificationMethod& method) -> void;

    /// Compute the thermodynamic properties of the phase.
    /// @param[in] props The evaluated thermodynamic properties of the phase.
    /// @param T The temperature of the phase (in K)
    /// @param P The pressure of the phase (in Pa)
    /// @param x The mole fractions of the species in the phase (in mol/mol)
    auto compute(CubicEOSProps& props, real T, real P, ArrayXrConstRef x) -> void;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
