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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/ReactionRateModel.hpp>

namespace Reaktoro {

/// The parameters in the reaction rate model of @cite{Palandri2004} for dissolution/precipitation kinetics of minerals.
struct ReactionRateModelParamsPalandriKharaka
{
    /// The parameters for a catalyser property in a mineral reaction rate mechanism of @cite{Palandri2004}.
    struct Catalyst
    {
        /// The chemical formula of the species that participates as a catalyst.
        String formula;

        /// The symbol of the species property that acts as a catalyser. The
        /// options are `a` for the activity of the species, and `P` for its
        /// partial pressure (in which case the species must be an existing gas
        /// in the system).
        String property = "a";

        /// The power of the property that affects the rate of mineral reaction.
        Param power = 0.0;
    };

    /// The parameters for a mineral reaction rate mechanism of @cite{Palandri2004}.
    struct Mechanism
    {
        /// The classifying name of the mineral reaction mechanism (e.g., `Acid`, `Neutral`, `Base`, `Carbonate`).
        String name;

        /// The kinetic rate constant of the mineral reaction at 298.15 K (in lg mol/(m2*s)).
        Param lgk;

        /// The Arrhenius activation energy of the mineral reaction (in kJ/mol).
        Param E;

        /// The empirical and dimensionless power parameter *p*.
        Param p = 1.0;

        /// The empirical and dimensionless power parameter *q*.
        Param q = 1.0;

        /// The catalysts of the mineral reaction.
        Vec<Catalyst> catalysts;
    };

    /// The name of the mineral (e.g., `Dolomite`).
    String mineral;

    /// The alternative names of the mineral (e.g., `Dolomite,ord`, `Dolomite,ordered`).
    Strings othernames;

    /// The reaction mechanisms considered in the mineral dissolution/precipitation rate model.
    Vec<Mechanism> mechanisms;
};

/// Return the reaction rate model of @cite{Palandri2004} for dissolution/precipitation kinetics of minerals.
/// The required model parameters will be fetched from the default parameters in `PalandriKharaka.yaml`.
auto ReactionRateModelPalandriKharaka() -> ReactionRateModelGenerator;

/// Return the reaction rate model of @cite{Palandri2004} for dissolution/precipitation kinetics of minerals.
/// The required model parameters will be fetched from the Params object `params`.
/// They must be available under a `PalandriKharaka` section inside a section `ReactionRateModelParams`.
/// @param params The object where mineral reaction rate parameters should be found.
auto ReactionRateModelPalandriKharaka(Params const& params) -> ReactionRateModelGenerator;

/// Return the reaction rate model of @cite{Palandri2004} for dissolution/precipitation kinetics of minerals.
auto ReactionRateModelPalandriKharaka(ReactionRateModelParamsPalandriKharaka const& params) -> ReactionRateModelGenerator;

/// Return the reaction rate model of @cite{Palandri2004} for dissolution/precipitation kinetics of minerals.
auto ReactionRateModelPalandriKharaka(Vec<ReactionRateModelParamsPalandriKharaka> const& paramsvec) -> ReactionRateModelGenerator;

} // namespace Reaktoro
