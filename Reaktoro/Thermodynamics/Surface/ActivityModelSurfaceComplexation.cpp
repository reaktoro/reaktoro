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
    const auto num_species = surface.species().size();

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

        // Calculate Davies only if the AqueousPhase has been already evaluated
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqstate = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Auxiliary constant references properties and variables
            const auto I = aqstate.Is;          // the stoichiometric ionic strength
            const auto sqrtI = sqrt(I);
            const auto ln10 = log(10);
            const auto Agamma = 0.5095;         // the Debye-Huckel parameter

            // Calculate the ln activity coefficient of the surface complexation species using the Davies activity model
            ln_g = ln10*(-Agamma*z*z*sqrtI/(1 + sqrtI) - 0.3*I);
        }

        // Add the correction introduced by the activity coefficients
        ln_a += ln_g;
    };
    return fn;
}

/// Return the SurfaceComplexationActivityModel object assuming the presence the Diffuse Double Layer (DDL) model.
auto activityModelSurfaceComplexationWithDDL(const SpeciesList& species, ActivityModelSurfaceComplexationParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

    // The number of surface complexation species in the current phase
    const auto num_species = surface.species().size();

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

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        ln_a = x.log();

        // Initialized the ln of activity coefficients of the surface complexation species
        ln_g = ArrayXr::Zero(num_species);

        // Auxiliary constant references properties
        real I;

        // Calculate the stoichiometric ionic strength if the DDL State has been already evaluated
        if (props.extra["DiffusiveLayerState"].has_value())
        {
            // Export ddl state via `extra` data member
            const auto& ddlstate = std::any_cast<AqueousMixtureState>(props.extra["DiffusiveLayerState"]);

            // Fetch the stoichiometric ionic strength
            I = ddlstate.Is;
        }
        // Otherwise, calculate the stoichiometric ionic strength if the Aqueous State has been already evaluated
        else if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqstate = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Fetch the stoichiometric ionic strength
            I = aqstate.Is;
        }

        // Evaluate surface complexation potential provided new ionic state
        surface_state.updatePotential(I);

        // Export the surface complexation and its state via the `extra` data member
        props.extra["ComplexationSurfaceState"] = surface_state;
        props.extra["ComplexationSurface"] = surface;

        // Auxiliary variables
        const auto sqrtI = sqrt(I);
        const auto ln10 = log(10);
        const auto Agamma = 0.5095; // the Debye-Huckel parameter

        // Calculate the ln activity coefficient of the surface complexation species using the Davies activity model
        ln_g = ln10*(-Agamma*z*z*sqrtI/(1 + sqrtI) - 0.3*I);

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
    const auto num_species = surface.species().size();

    // The equivalent numbers (absolute values of charges) of the surface complexation species
    ArrayXd ze = surface.equivalentsNumbers();

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

            // Auxiliary constant references properties and variables
            const auto I = aqstate.Is;            // the stoichiometric ionic strength
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

/// Return the Diffuse Double Layer (DDL) activity model based on the Dzombak and Morel (1990) model
auto activityModelDDL(const SpeciesList& species, ActivityModelDDLParams params) -> ActivityModel
{
    // Create the aqueous ddl_mixture
    AqueousMixture ddl_mixture(species);

    // The number of all species and surface complexation species in the current exchange phase only
    const auto num_species = species.size();

    // The ddl_state of the aqueous ddl_mixture
    AqueousMixtureState ddl_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the ddl_state of the aqueous ddl_mixture
        ddl_state = ddl_mixture.state(T, P, x);

        // Export the surface complexation and its ddl_state via the `extra` data member
        props.extra["DiffusiveLayerState"] = ddl_state;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        //  ln_a = log(params.enr_dll) + log(ddl_state.m);

        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export surface complexation ddl_state via `extra` data member
            const auto& aq_state = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Calculate ln(a) of the DDL layer species, according to the formula:
            // cD = c*enr,
            // where c   is the concentration of the species in the aqueous solution and
            //       enr is the enrichment factor.
            ln_a = log(params.enr_dll) + log(aq_state.m);

            const auto& aq_mix = std::any_cast<AqueousMixture>(props.extra["AqueousMixture"]);
            const auto& z_aq = aq_mix.charges();
            const auto& species_aq = aq_mix.species();
            std::cout << "aq_state   = " << aq_state.ms.transpose() << std::endl;
            std::cout << "z_aq       = " << z_aq.transpose() << std::endl;
            for(auto s : species_aq) std::cout << s.name() << " ";
            std::cout << std::endl;
            const auto Zaq = (z_aq * aq_state.m).sum();
            std::cout << "Zaq = " << Zaq << std::endl;
            getchar();

        }
        if (props.extra["ComplexationSurfaceState"].has_value())
        {
            // Export surface complexation ddl_state via `extra` data member
            auto surface_state = std::any_cast<ComplexationSurfaceState>(props.extra["ComplexationSurfaceState"]);

            // Auxiliary constants
            const auto F = faradayConstant;
            const auto R = universalGasConstant;

            // Auxiliary constant references properties of the surface
            const auto& sigma = surface_state.sigma; // the electrostatic surface potential
            const auto& Z = surface_state.Z;        // the charges of the surface species

            // Auxiliary constant references properties of the DDL
            const auto I = ddl_state.Ie;
            surface_state.updatePotential(I);

            // Update activity coefficient using the coulombic correction factor
            ln_g = -Z*F*surface_state.psi/(R*T) * ArrayXr::Ones(num_species);
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

auto ActivityModelSurfaceComplexationWithDDL(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationWithDDL(surface_species, params);
    };
}

auto ActivityModelSurfaceComplexationGainesThomas(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationGainesThomas(surface_species, params);
    };
}

auto ActivityModelDDL() -> ActivityModelGenerator
{
    ActivityModelDDLParams params;

    return [=](const SpeciesList& species)
    {
        return detail::activityModelDDL(species, params);
    };
}

auto ActivityModelDDL(ActivityModelDDLParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return detail::activityModelDDL(species, params);
    };
}

} // namespace Reaktoro
