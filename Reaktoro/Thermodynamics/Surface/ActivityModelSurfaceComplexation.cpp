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

#include "ActivityModelSurfaceComplexation.hpp"

// Reaktoro includes
#include <Reaktoro/Singletons/Elements.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

using std::sqrt;
using std::log;

namespace detail {

/// Return the SurfaceComplexationActivityModel object assuming no electrostatic effects.
auto activityModelSurfaceComplexationNoDDL(const SpeciesList& species, ActivityModelSurfaceComplexationParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

    // The number of surface complexation species in the current phase
    const auto num_species = species.size();

    // The charges of the surface complexation species
    ArrayXd z = surface.charges();

    // The state of the complexation surface
    ComplexationSurfaceState surface_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the surface complexation
        surface_state = surface.state(T, P, x);

        // Export the surface complexation and its state via the `extra` data member
        props.extra["ComplexationSurfaceState"] = surface_state;
        props.extra["ComplexationSurface"] = surface;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        ln_a = x.log();

        // Initialized the ln of activity coefficients of the surface complexation species
        ln_g = ArrayXr::Zero(num_species);

        // Calculate Davies and Debye-Huckel parameters only if the AqueousPhase has been already evaluated
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqstate = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Auxiliary constant references properties
            const auto& I = aqstate.Is;            // the stoichiometric ionic strength
            const auto& rho = aqstate.rho/1000;    // the density of water (in g/cm3)
            const auto& epsilon = aqstate.epsilon; // the dielectric constant of water

            // Auxiliary variables
            const auto sqrtI = sqrt(I);
            const auto ln10 = log(10);
            const auto Agamma = 0.5095; // the Debye-Huckel parameter

            // Calculate the ln activity coefficient of the surface complexation species using the Davies activity model
            ln_g = ln10*(-Agamma*z*z*sqrtI/(1 + sqrtI) - 0.3*I);
        }

        // Add the correction introduced by the activity coefficients
        ln_a += ln_g;
    };
    return fn;
}

/// Return the SurfaceComplexationActivityModel object assuming no electrostatic effects and the Gaines-Thomas convention
/// for the activities calculation.
auto activityModelSurfaceComplexationGainesThomas(const SpeciesList& species, ActivityModelSurfaceComplexationParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

    // The number of surface complexation species in the current phase
    const auto num_species = species.size();

    // The equivalent numbers (absolute values of charges) of the surface complexation species
    ArrayXd ze = abs(surface.charges());

    // The state of the complexation surface
    ComplexationSurfaceState surface_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the surface complexation
        surface_state = surface.state(T, P, x);

        // Export the surface complexation and its state via the `extra` data member
        props.extra["ComplexationSurfaceState"] = surface_state;
        props.extra["ComplexationSurface"] = surface;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Calculate ln of activities of surfaces species as the ln of equivalent fractions
        ln_a = (x*ze/(x*ze).sum()).log();

        // Initialized the ln of activity coefficients of the surface complexation species
        ln_g = ArrayXr::Zero(num_species);

        // Calculate Davies and Debye-Huckel parameters only if the AqueousPhase has been already evaluated
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqstate = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Auxiliary constant references properties
            const auto& I = aqstate.Is;            // the stoichiometric ionic strength
            const auto& rho = aqstate.rho/1000;    // the density of water (in g/cm3)
            const auto& epsilon = aqstate.epsilon; // the dielectric constant of water

            // Auxiliary variables
            const auto sqrtI = sqrt(I);
            const auto ln10 = log(10);
            const auto Agamma = 0.5095; // the Debye-Huckel parameter

            // Calculate the ln activity coefficient of the surface complexation species using the Davies activity model
            ln_g = ln10*(-Agamma*ze*ze*sqrtI/(1 + sqrtI) - 0.3*I);
        }

        // Add the correction introduced by the activity coefficients
        ln_a += ln_g;
    };
    return fn;
}

/// Return the SurfaceComplexationActivityModel object based on the Diffuse Double Layer (DDL) model.
auto activityModelSurfaceComplexationDDL(const SpeciesList& species, ActivityModelSurfaceComplexationDDLParams params) -> ActivityModel
{
    // Create the aqueous mixture
    AqueousMixture mixture(species);

    // The number of all species and surface complexation species in the current exchange phase only
    const auto num_species = species.size();

    // The state of the aqueous mixture
    AqueousMixtureState state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the aqueous mixture
        state = mixture.state(T, P, x);

        // Export the surface complexation and its state via the `extra` data member
        props.extra["DiffusiveLayerState"] = state;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Calculate ln of activities of aqueous species in the diffusive layer
        // as the ln of molar fractions multiplied by the enrichment factor
        ln_a = params.enrichment_factor_dll*x.log();

        // Initialized the ln of activity coefficients of the surface complexation species on the surface
        ln_g = ArrayXr::Zero(num_species);

        // Calculate Davies and Debye--Huckel parameters only if the AqueousPhase has been already evaluated
        if (props.extra["SurfaceComplexationState"].has_value())
        {
            // Export surface complexation state via `extra` data member
            const auto& surfcomplex_state = std::any_cast<ComplexationSurfaceState>(props.extra["ComplexationSurfaceState"]);

            // Auxiliary constants
            const auto F = faradayConstant;
            const auto R = universalGasConstant;

            // Auxiliary constant references properties
            const auto& sigma = surfcomplex_state.sigma; // the electrostatic surface potential
            const auto& z = surfcomplex_state.z;        // the charges of the surface species
            const auto I = state.Ie;

            // Using formula sigma = 0.1174*I^0.5*sinh(F*psi/R/T/2) and arcsinh(y) = ln(y+(y^2+1)^1⁄2)
            const auto y = sigma/(0.1174*sqrt(I));
            const auto psi = 2*R*T/F*log(y + sqrt(1 + y*y));

            // Calculate the lng according to the formula g = exp(-z*F*psi/R/T)
            ln_g = -z*F*psi/(R*T);
        }

        // Add the correction introduced by the activity coefficients
        ln_a += ln_g;
    };
    return fn;
}

} // namespace detail

auto ActivityModelSurfaceComplexationNoDDL(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationNoDDL(surface_species, params);
    };
}

auto ActivityModelSurfaceComplexationGainesThomas(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationGainesThomas(surface_species, params);
    };
}

auto ActivityModelSurfaceComplexationDDL() -> ActivityModelGenerator
{
    ActivityModelSurfaceComplexationDDLParams params;

    return [=](const SpeciesList& species)
    {
        return detail::activityModelSurfaceComplexationDDL(species,  params);
    };
}

auto ActivityModelSurfaceComplexationDDL(ActivityModelSurfaceComplexationDDLParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return detail::activityModelSurfaceComplexationDDL(species,  params);
    };
}

} // namespace Reaktoro
