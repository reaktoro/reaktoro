// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/Params.hpp>

namespace Reaktoro {

/// The parameters in the Extended UNIQUAC activity model for aqueous electrolyte solutions.
struct ActivityModelParamsExtendedUNIQUAC
{
    /// The available temperature-pressure correction models for the Extended UNIQUAC interaction parameters.
    enum class CorrectionModel
    {
        Constant,  ///< Correction model based on formula \eq{\mathrm{Param}(T,P)=c_{0}}.
        Linear,    ///< Correction model based on formula \eq{\mathrm{Param}(T,P)=c_{0}+c_{1}T}.
        Quadratic, ///< Correction model based on formula \eq{\mathrm{Param}(T,P)=c_{0}+c_{1}T+c_{1}T^2}.
    };

    /// The attributes of an interaction parameter in the Extended UNIQUAC activity model.
    struct InteractionParamAttribs
    {
        /// The chemical formulas of the species associated to this species interaction parameter sorted in descending order of charge.
        Vec<ChemicalFormula> formulas;

        /// The model used for temperature-pressure correction of this species interaction parameter (options).
        CorrectionModel model = CorrectionModel::Constant;

        /// The parameters for the temperature-pressure correction model of this species interaction parameter.
        Vec<Param> parameters;
    };

    Vec<InteractionParamAttribs> beta0;  ///< The parameters \eq{\beta^{(0)}_{ij}(T, P)} in the Extended UNIQUAC model for cation-anion interactions.
    Vec<InteractionParamAttribs> beta1;  ///< The parameters \eq{\beta^{(1)}_{ij}(T, P)} in the Extended UNIQUAC model for cation-anion interactions.
    Vec<InteractionParamAttribs> beta2;  ///< The parameters \eq{\beta^{(2)}_{ij}(T, P)} in the Extended UNIQUAC model for cation-anion interactions.
    Vec<InteractionParamAttribs> Cphi;   ///< The parameters \eq{C^{\phi}_{ij}(T, P)} in the Extended UNIQUAC model for cation-anion interactions.
    Vec<InteractionParamAttribs> theta;  ///< The parameters \eq{\theta_{ij}(T, P)} in the Extended UNIQUAC model for cation-cation and anion-anion interactions.
    Vec<InteractionParamAttribs> lambda; ///< The parameters \eq{\lambda_{ij}(T, P)} in the Extended UNIQUAC model for neutral-cation and neutral-anion interactions.
    Vec<InteractionParamAttribs> zeta;   ///< The parameters \eq{\zeta_{ijk}(T, P)} in the Extended UNIQUAC model for neutral-cation-anion interactions.
    Vec<InteractionParamAttribs> psi;    ///< The parameters \eq{\psi_{ijk}(T, P)} in the Extended UNIQUAC model for cation-cation-anion and anion-anion-cation interactions.
    Vec<InteractionParamAttribs> mu;     ///< The parameters \eq{\mu_{ijk}(T, P)} in the Extended UNIQUAC model for neutral-neutral-neutral, neutral-neutral-cation, and neutral-neutral-anion interactions.
    Vec<InteractionParamAttribs> eta;    ///< The parameters \eq{\eta_{ijk}(T, P)} in the Extended UNIQUAC model for neutral-cation-cation and neutral-anion-anion interactions.
    Vec<InteractionParamAttribs> alpha1; ///< The parameters \eq{\alpha_1_{ij}} associated to the parameters \eq{\beta^{(1)}_{ij}}.
    Vec<InteractionParamAttribs> alpha2; ///< The parameters \eq{\alpha_2_{ij}} associated to the parameters \eq{\beta^{(2)}_{ij}}.
};

/// Return the activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelExtendedUNIQUAC() -> ActivityModelGenerator;

/// Return the activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelExtendedUNIQUAC(ActivityModelParamsExtendedUNIQUAC const& params) -> ActivityModelGenerator;

/// Return the activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelExtendedUNIQUAC(Params const& params) -> ActivityModelGenerator;

} // namespace Reaktoro
