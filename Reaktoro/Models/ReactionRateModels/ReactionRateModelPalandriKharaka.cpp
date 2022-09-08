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

#include "ReactionRateModelPalandriKharaka.hpp"

// Reaktoro includes
// #include <Reaktoro/Common/Constants.hpp>
// #include <Reaktoro/Serialization/Models.YAML.hpp>

namespace Reaktoro {

auto ReactionRateModelPalandriKharaka(ReactionRateModelParamsPalandriKharaka const& params) -> ReactionRateModelGenerator
{
    ReactionRateModelGenerator model = [](Reaction const& reaction, PhaseList const& phases)
    {
        // const auto iH2O = species.indexWithFormula("H2O");
        // const auto MH2O = waterMolarMass;

        ReactionRateModel fn = [=](ChemicalProps const& props) -> Rate
        {
            // const auto x = args.x;
            // const auto m = x/(MH2O * args.x[iH2O]); // molalities

            // // Set the state of matter of the phase
            // props.som = StateOfMatter::Liquid;

            // props = 0.0;
            // props.ln_a = m.log();
            // props.ln_a[iH2O] = log(x[iH2O]);
            return {};
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
