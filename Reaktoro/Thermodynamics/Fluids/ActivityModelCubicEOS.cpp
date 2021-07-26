// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/StateOfMatter.hpp>
#include <Reaktoro/Singletons/CriticalProps.hpp>
#include <Reaktoro/Thermodynamics/Fluids/CubicEOS.hpp>

namespace Reaktoro {

using std::log;

auto activityPropsFnCubicEOS(const SpeciesList& species, ActivityModelCubicEOSParams params, CubicEOSModel type) -> ActivityPropsFn
{
    // The number of gases
    const auto nspecies = species.size();

    // Get the the critical temperatures, pressures and acentric factors of the gases
    ArrayXr Tcr(nspecies), Pcr(nspecies), omega(nspecies);
    for(auto i = 0; i < nspecies; ++i)
    {
        const auto crprops = CriticalProps::get({
            species[i].substance(),
            species[i].formula(),
            species[i].name()
        });
        error(!crprops.has_value(), "Cannot create any cubic equation of state model "
            "(e.g. Peng-Robinson, Soave-Redlich-Kwong, etc.) without "
            "critical properties for the species with name ", species[i].name(), ". "
            "In order to fix this error, use CriticalProps::append to register the "
            "critical properties of this substance.");
        Tcr[i] = crprops->temperature();
        Pcr[i] = crprops->pressure();
        omega[i] = crprops->acentricFactor();
    }

    const auto aggregatestate = species[0].aggregateState();

    error(aggregatestate != AggregateState::Gas && aggregatestate != AggregateState::Liquid,
        "Cannot create a cubic equation of state model if the species "
        "in the phase have aggregate state ", aggregatestate, ". "
        "Only Gas or Liquid AggregateState values are permitted.");

    // Initialize the CubicEOS instance
    CubicEOS eos({nspecies, Tcr, Pcr, omega});

    if( aggregatestate == AggregateState::Gas )
        eos.setFluidType(CubicEOSFluidType::Vapor);
    else eos.setFluidType(CubicEOSFluidType::Liquid);

    eos.setModel(type);
    eos.setInteractionParamsFunction(params.interaction_params_fn);
    eos.setStablePhaseIdentificationMethod(params.phase_identification_method);

    /// The thermodynamic properties calculated with CubicEOS
    CubicEOSProps res;
    res.ln_phi.resize(nspecies);

    // Define the activity model function of the gaseous phase
    ActivityPropsFn fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x, extra] = args;

        const auto Pbar = P * 1.0e-5; // convert from Pa to bar

        eos.compute(res, T, P, x);

        props.Vex  = res.Vres;
        props.VexT = res.VresT;
        props.VexP = res.VresP;
        props.Gex  = res.Gres;
        props.Hex  = res.Hres;
        props.Cpex = res.Cpres;
        props.Cvex = res.Cvres;
        props.ln_g = res.ln_phi;
        props.ln_a = res.ln_phi + log(x) + log(Pbar);
    };

    return fn;
}

auto ActivityModelCubicEOS(ActivityModelCubicEOSParams params, CubicEOSModel type) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return activityPropsFnCubicEOS(species, params, type);
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
