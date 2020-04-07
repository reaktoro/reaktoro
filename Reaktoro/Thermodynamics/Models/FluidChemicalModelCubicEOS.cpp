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

auto activityModelCubicEOS(const GeneralMixture& mixture, ActivityModelOptionsCubicEOS options) -> ActivityModelFn
{
    // The number of gases in the mixture
    const auto nspecies = mixture.numSpecies();

    // Get the the critical temperatures, pressures and acentric factors of the gases
    ArrayXr Tcr(nspecies), Pcr(nspecies), omega(nspecies);
    for(auto i = 0; i < nspecies; ++i)
    {
        const auto& species = mixture.species()[i];
        const auto crprops = species.criticalProps();
        error(!crprops.has_value(), "Cannot create any cubic equation of state model (e.g. Peng-Robinson, Soave-Redlich-Kwong, etc.) without "
            "critical properties for the species with symbol ", species.symbol(), "and with substance name", species.name(), ". "
            "In order to fix this error, use CriticalProps::append to register the critical properties of this substance.");
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
    CubicEOSProps props;
    props.ln_phi.resize(nspecies);

    // Define the chemical model function of the gaseous phase
    ActivityModelFn model = [=](ActivityProps& res, real T, real P, ArrayXrConstRef x) mutable
    {
        // Evaluate the CubicEOS
        eos.compute(props, T, P, x);

        // The ln of mole fractions
        const auto ln_x = log(x);

        // The ln of pressure in bar units
        const auto ln_Pbar = log(1e-5 * P);

        // Fill the chemical properties of the fluid phase
        res.Vex  = props.Vres;
        res.VexT = props.VresT;
        res.VexP = props.VresP;
        res.Gex  = props.Gres;
        res.Hex  = props.Hres;
        res.Cpex = props.Cpres;
        res.ln_g = props.ln_phi;
        res.ln_a = props.ln_phi + ln_x + ln_Pbar;
    };

    return model;
}

} // namespace Reaktoro
