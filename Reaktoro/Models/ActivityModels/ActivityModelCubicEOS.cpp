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
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Models/ActivityModels/Support/CubicEOS.hpp>
#include <Reaktoro/Singletons/CriticalProps.hpp>

namespace Reaktoro {

using std::log;

/// Convert the list of Species objects into a list of CubicEOS::Substance objects for CubicEOS::Equation.
auto collectSubstancesInFluidPhase(SpeciesList const& specieslist) -> Vec<CubicEOS::Substance>
{
    Vec<CubicEOS::Substance> substances;
    substances.reserve(specieslist.size());

    for(auto species : specieslist)
    {
        const auto crprops = CriticalProps::get({
            species.substance(),
            species.formula(),
            species.name()
        });

        errorifnot(crprops.has_value(),
            "Could not find critical properties for substance ", species.formula().str(), " "
            "while creating a cubic equation of state model (e.g. Peng-Robinson, Soave-Redlich-Kwong, etc.). ",
            "The following are ways to fix this error:\n"
            "  - use CriticalProps::append(crprops) to register the critical properties of this substance in the CriticalProps database,\n"
            "  - use CriticalProps::setMissingAs(\"He\") to consider all missing substances in the CriticalProps database as if they were He.");


        auto const& formula = species.formula();
        auto const& Tcr = crprops->temperature();
        auto const& Pcr = crprops->pressure();
        auto const& omega = crprops->acentricFactor();

        substances.push_back(CubicEOS::Substance{formula, Tcr, Pcr, omega});
    }

    return substances;
}

auto activityModelCubicEOS(SpeciesList const& specieslist, CubicEOS::EquationModel const& eqmodel, CubicEOS::BipModel const& bipmodel) -> ActivityModel
{
    CubicEOS::EquationSpecs eqspecs;
    eqspecs.substances = collectSubstancesInFluidPhase(specieslist);
    eqspecs.eqmodel = eqmodel;
    eqspecs.bipmodel = bipmodel;

    CubicEOS::Equation equation(eqspecs);

    CubicEOS::Props cprops;

    // Define the activity model function of the fluid phase
    ActivityModel model = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& [T, P, x] = args;

        const auto Pbar = P * 1.0e-5; // convert from Pa to bar

        equation.compute(cprops, T, P, x);

        props.Vx   = cprops.V;
        props.VxT  = cprops.VT;
        props.VxP  = cprops.VP;
        props.Gx   = cprops.Gres;
        props.Hx   = cprops.Hres;
        props.Cpx  = cprops.Cpres;
        props.ln_g = cprops.ln_phi;
        props.ln_a = cprops.ln_phi + log(x) + log(Pbar);
        props.som  = cprops.som;
    };

    return model;
}

auto CubicBipModelPhreeqc() -> CubicBipModelGenerator
{
    return [](SpeciesList const& specieslist) -> CubicEOS::BipModel
    {
        Strings substances = vectorize(specieslist, RKT_LAMBDA(x, x.formula().str()));
        return CubicEOS::BipModelPhreeqc(substances, CubicEOS::BipModelParamsPhreeqc());
    };
}

auto CubicBipModelSoreideWhitson() -> CubicBipModelGenerator
{
    return [](SpeciesList const& specieslist) -> CubicEOS::BipModel
    {
        Strings substances = vectorize(specieslist, RKT_LAMBDA(x, x.formula().str()));
        return CubicEOS::BipModelSoreideWhitson(substances, CubicEOS::BipModelParamsSoreideWhitson());
    };
}

auto ActivityModelCubicEOS(CubicBipModelGenerator cbipmodel, CubicEOS::EquationModel const& eqmodel) -> ActivityModelGenerator
{
    return [=](SpeciesList const& specieslist) -> ActivityModel
    {
        CubicEOS::BipModel bipmodel = cbipmodel ? cbipmodel(specieslist) : CubicEOS::BipModel{};
        return activityModelCubicEOS(specieslist, eqmodel, bipmodel);
    };
}

auto ActivityModelVanDerWaals(CubicBipModelGenerator cbipmodel) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(cbipmodel, CubicEOS::EquationModelVanDerWaals());
}

auto ActivityModelRedlichKwong(CubicBipModelGenerator cbipmodel) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(cbipmodel, CubicEOS::EquationModelRedlichKwong());
}

auto ActivityModelSoaveRedlichKwong(CubicBipModelGenerator cbipmodel) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(cbipmodel, CubicEOS::EquationModelSoaveRedlichKwong());
}

auto ActivityModelPengRobinson(CubicBipModelGenerator cbipmodel) -> ActivityModelGenerator
{
    return ActivityModelPengRobinson78(cbipmodel);
}

auto ActivityModelPengRobinson76(CubicBipModelGenerator cbipmodel) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(cbipmodel, CubicEOS::EquationModelPengRobinson76());
}

auto ActivityModelPengRobinson78(CubicBipModelGenerator cbipmodel) -> ActivityModelGenerator
{
    return ActivityModelCubicEOS(cbipmodel, CubicEOS::EquationModelPengRobinson78());
}

auto ActivityModelPengRobinsonPhreeqc() -> ActivityModelGenerator
{
    return ActivityModelPengRobinson76(CubicBipModelPhreeqc());
}

auto ActivityModelPengRobinsonSoreideWhitson() -> ActivityModelGenerator
{
    return ActivityModelPengRobinson78(CubicBipModelSoreideWhitson());
}

// METHODS BELOW ARE DEPRECATED AND WILL BE REMOVED IN FUTURE RELEASES

auto ActivityModelPengRobinsonPHREEQC() -> ActivityModelGenerator
{
    return ActivityModelPengRobinsonPhreeqc();
}

auto CubicBipModelPHREEQC() -> CubicBipModelGenerator
{
    return CubicBipModelPhreeqc();
}

} // namespace Reaktoro
