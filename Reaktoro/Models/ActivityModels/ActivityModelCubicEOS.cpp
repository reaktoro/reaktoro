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

#include "ActivityModelCubicEOS.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Singletons/CriticalProps.hpp>
#include <Reaktoro/Models/ActivityModels/Support/CubicEOS.hpp>

namespace Reaktoro {

using std::log;

auto activityModelCubicEOS(const SpeciesList& species, ActivityModelCubicEOSParams params, CubicEOSModel type) -> ActivityModel
{
    // The number of gases
    const auto nspecies = species.size();

    // Get the critical temperatures, pressures and acentric factors of the gases
    ArrayXr Tcr(nspecies), Pcr(nspecies), omega(nspecies);
    for(auto i = 0; i < nspecies; ++i)
    {
        const auto crprops = CriticalProps::get({
            species[i].substance(),
            species[i].formula(),
            species[i].name()
        });

        errorif(crprops.has_value() == false,
            "Could not find critical properties for substance ", species[i].formula().str(), " "
            "while creating a cubic equation of state model (e.g. Peng-Robinson, Soave-Redlich-Kwong, etc.). ",
            "The following are ways to fix this error:\n"
            "  - use CriticalProps::append(crprops) to register the critical properties of this substance in the CriticalProps database,\n"
            "  - use CriticalProps::setMissingAs(\"He\") to consider all missing substances in the CriticalProps database as if they were He.");

        Tcr[i] = crprops->temperature();
        Pcr[i] = crprops->pressure();
        omega[i] = crprops->acentricFactor();
    }

    // Initialize the CubicEOS instance
    CubicEOS eos({nspecies, Tcr, Pcr, omega});

    eos.setModel(type);
    eos.setInteractionParamsFunction(params.interaction_params_fn);

    /// The thermodynamic properties calculated with CubicEOS
    CubicEOSProps res;
    res.ln_phi.resize(nspecies);

    // Define the activity model function of the gaseous phase
    ActivityModel model = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        const auto Pbar = P * 1.0e-5; // convert from Pa to bar

        eos.compute(res, T, P, x);

        props.Vx   = res.V;
        props.VxT  = res.VT;
        props.VxP  = res.VP;
        props.Gx   = res.Gres;
        props.Hx   = res.Hres;
        props.Cpx  = res.Cpres;
        props.ln_g = res.ln_phi;
        props.ln_a = res.ln_phi + log(x) + log(Pbar);
        props.som  = res.som;
    };

    return model;
}

auto ActivityModelCubicEOS(ActivityModelCubicEOSParams params, CubicEOSModel type) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return activityModelCubicEOS(species, params, type);
    };
}

auto ActivityModelVanDerWaals(ActivityModelCubicEOSParams params) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(params, CubicEOSModel::VanDerWaals);
}

auto ActivityModelRedlichKwong(ActivityModelCubicEOSParams params) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(params, CubicEOSModel::RedlichKwong);
}

auto ActivityModelSoaveRedlichKwong(ActivityModelCubicEOSParams params) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(params, CubicEOSModel::SoaveRedlichKwong);
}

auto ActivityModelPengRobinson(ActivityModelCubicEOSParams params) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(params, CubicEOSModel::PengRobinson);
}

} // namespace Reaktoro
