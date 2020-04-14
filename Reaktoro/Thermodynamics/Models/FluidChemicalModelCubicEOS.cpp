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

#include "FluidChemicalModelCubicEOS.hpp"

// C++ includes
#include <map>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>
#include <Reaktoro/Singletons/CriticalProps.hpp>
#include <Reaktoro/Thermodynamics/EOS/CubicEOS.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GeneralMixture.hpp>

namespace Reaktoro {

using std::log;

auto activityModelCubicEOS(const GeneralMixture& mixture, ActivityModelOptionsCubicEOS options) -> ActivityPropsFn
{
    // The number of gases in the mixture
    const auto nspecies = mixture.numSpecies();

    // Get the the critical temperatures, pressures and acentric factors of the gases
    ArrayXr Tcr(nspecies), Pcr(nspecies), omega(nspecies);
    for(auto i = 0; i < nspecies; ++i)
    {
        const auto& species = mixture.species()[i];
        const auto crprops = species.criticalProps();
        error(!crprops.has_value(), "Cannot create any cubic equation of state model "
            "(e.g. Peng-Robinson, Soave-Redlich-Kwong, etc.) without "
            "critical properties for the species with name ", species.name(), ". "
            "In order to fix this error, use CriticalProps::append to register the "
            "critical properties of this substance.");
        Tcr[i] = crprops->temperature();
        Pcr[i] = crprops->pressure();
        omega[i] = crprops->acentricFactor();
    }

    // Initialize the CubicEOS instance
    CubicEOS eos({nspecies, Tcr, Pcr, omega});
    eos.setFluidType(options.fluidtype);
    eos.setModel(options.model);
    eos.setInteractionParamsFunction(options.interaction_params_fn);
    eos.setStablePhaseIdentificationMethod(options.phase_identification_method);

    /// The thermodynamic properties calculated with CubicEOS
    CubicEOSProps res;
    res.ln_phi.resize(nspecies);

    // Define the activity model function of the gaseous phase
    ActivityPropsFn fn = [=](ActivityProps props, real T, real P, ArrayXrConstRef x) mutable
    {
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
        props.ln_a = res.ln_phi + log(x) + log(P);
    };

    return fn;
}

} // namespace Reaktoro
