// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

/// The parameters in the Pitzer activity model for aqueous electrolyte solutions.
struct ActivityModelParamsPitzer
{
    /// The available temperature-pressure correction models for the Pitzer interaction parameters.
    enum class CorrectionModel
    {
        Constant,           ///< Correction model based on formula \eq{\mathrm{Param}(T,P)=c_{0}}.
        Phreeqc,            ///< Correction model based on formula \eq{\mathrm{Param}(T,P)=c_{0}+c_{1}\left(\dfrac{1}{T}-\dfrac{1}{T_{r}}\right)+c_{2}\ln\left(\dfrac{T}{T_{r}}\right)+c_{3}(T-T_{r})+c_{4}(T^{2}-T_{r}^{2})+c_{5}\left(\dfrac{1}{T^{2}}-\dfrac{1}{T_{r}^{2}}\right)} from PHREEQC v3 where \eq{T_{r}=298.15\;\text{K}}.
        HeMorse1993,        ///< Correction model based on formula \eq{\mathrm{Param}(T,P)=c_{0}+c_{1}T+c_{2}T^{2}+c_{3}/T+c_{4}\ln(T)} from He and Morse (1993) (doi: 10.1016/0016-7037(93)90137-L).
        Dai2013,            ///< Correction model based on formula \eq{\mathrm{Param}(T,P)=c_{0}+\dfrac{c_{1}}{T}+c_{2}T+c_{3}\ln T+c_{4}T^{2}+\left[c_{5}+\dfrac{c_{6}}{T}+c_{7}T+c_{8}\ln T+c_{9}T^{2}\right]P+\left[c_{10}+\dfrac{c_{11}}{T}+c_{12}T+c_{13}\ln T+c_{14}T^{2}\right]P^{2}+\left[c_{15}+\dfrac{c_{16}}{T}+c_{17}T+c_{18}\ln T+c_{19}T^{2}\right]P^{3}} from Dai et al. (2013) (doi: 10.2118/164045-ms).
        Dai2014,            ///< Correction model based on formula \eq{\mathrm{Param}(T,P)=c_{0}+\dfrac{c_{1}}{T}+c_{2}T+c_{3}\ln T+\left[c_{4}+\dfrac{c_{5}}{T}+c_{6}T+c_{7}\ln T\right]P+\left[c_{8}+\dfrac{c_{9}}{T}+c_{10}T+c_{11}\ln T\right]P^{2}} from Dai et al. (2014) (doi: 10.2118/169786-ms).
        ChristovMoller2004, ///< Correction model based on formula \eq{} from Christov and Møller (2004) (see Table 5 in 10.1016/j.chemgeo.2007.07.023).
        Holmes1987,         ///< Correction model based on formula \eq{} from Holmes et al. (1987) (see Table 6 in 10.1016/j.chemgeo.2007.07.023).
        Pitzer1984,         ///< Correction model based on formula \eq{} from Pitzer et al. (1984) (see Table 7 in 10.1016/j.chemgeo.2007.07.023).
        PalabanPitzer1987,  ///< Correction model based on formula \eq{} from Palaban and Pitzer (1987) (see Table 8 in 10.1016/j.chemgeo.2007.07.023).
        Polya2001,          ///< Correction model based on formula \eq{} from Polya et al. (2001) (see Table 10 in 10.1016/j.chemgeo.2007.07.023).
        LiDuan2007,         ///< Correction model based on formula \eq{} from Li and Duan (2007) (see Table 12 in 10.1016/j.chemgeo.2007.07.023).
    };

    /// The attributes of an interaction parameter in the Pitzer activity model.
    struct InteractionParamAttribs
    {
        /// The chemical formulas of the species associated to this species interaction parameter sorted in descending order of charge.
        Vec<ChemicalFormula> formulas;

        /// The model used for temperature-pressure correction of this species interaction parameter (options).
        CorrectionModel model = CorrectionModel::Constant;

        /// The parameters for the temperature-pressure correction model of this species interaction parameter.
        Vec<real> parameters;
    };

    Vec<InteractionParamAttribs> beta0;  ///< The parameters \eq{\beta^{(0)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<InteractionParamAttribs> beta1;  ///< The parameters \eq{\beta^{(1)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<InteractionParamAttribs> beta2;  ///< The parameters \eq{\beta^{(2)}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<InteractionParamAttribs> Cphi;   ///< The parameters \eq{C^{\phi}_{ij}(T, P)} in the Pitzer model for cation-anion interactions.
    Vec<InteractionParamAttribs> theta;  ///< The parameters \eq{\theta_{ij}(T, P)} in the Pitzer model for cation-cation and anion-anion interactions.
    Vec<InteractionParamAttribs> lambda; ///< The parameters \eq{\lambda_{ij}(T, P)} in the Pitzer model for neutral-cation and neutral-anion interactions.
    Vec<InteractionParamAttribs> zeta;   ///< The parameters \eq{\zeta_{ijk}(T, P)} in the Pitzer model for neutral-cation-anion interactions.
    Vec<InteractionParamAttribs> psi;    ///< The parameters \eq{\psi_{ijk}(T, P)} in the Pitzer model for cation-cation-anion and anion-anion-cation interactions.
    Vec<InteractionParamAttribs> mu;     ///< The parameters \eq{\mu_{ijk}(T, P)} in the Pitzer model for neutral-neutral-neutral, neutral-neutral-cation, and neutral-neutral-anion interactions.
    Vec<InteractionParamAttribs> eta;    ///< The parameters \eq{\eta_{ijk}(T, P)} in the Pitzer model for neutral-cation-cation and neutral-anion-anion interactions.
    Vec<InteractionParamAttribs> alpha1; ///< The parameters \eq{\alpha_1_{ij}} associated to the parameters \eq{\beta^{(1)}_{ij}}.
    Vec<InteractionParamAttribs> alpha2; ///< The parameters \eq{\alpha_2_{ij}} associated to the parameters \eq{\beta^{(2)}_{ij}}.
};

/// Return the Pitzer activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelPitzer() -> ActivityModelGenerator;

/// Return the Pitzer activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelPitzer(ActivityModelParamsPitzer const& params) -> ActivityModelGenerator;

/// Return the Pitzer activity model for aqueous electrolyte phases based on PHREEQC's implementation.
auto ActivityModelPitzer(Params const& params) -> ActivityModelGenerator;

} // namespace Reaktoro
