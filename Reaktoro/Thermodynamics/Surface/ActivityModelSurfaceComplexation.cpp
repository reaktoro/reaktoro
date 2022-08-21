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
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurfaceSite.hpp>
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

using std::sqrt;
using std::log;

namespace detail {

/// Return the SurfaceComplexationActivityModel object assuming no electrostatic effects and considering every site as
/// a separate phase.
auto activityModelSurfaceComplexationSiteNoDDL(const SpeciesList& species, ActivityModelSurfaceComplexationSiteParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

    // Create the complexation surface site
    ComplexationSurfaceSite site = surface.sites()[params.site_tag];

    // The state of the complexation surface
    ComplexationSurfaceSiteState site_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the surface complexation site
        site_state = site.state(T, P, x);

        // Export the surface complexation and its state via the `extra` data member
        props.extra["ComplexationSurfaceSiteState" + site.name()] = site_state;
        props.extra["ComplexationSurfaceSite" + site.name()] = site;

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        props.ln_a = x.log();
    };

    return fn;
}

/// Return the SurfaceComplexationActivityModel object assuming no electrostatic effects.
auto activityModelSurfaceComplexationNoDDL(const SpeciesList& species, ActivityModelSurfaceComplexationParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

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

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        props.ln_a = x.log();
    };

    return fn;
}

/// Return the SurfaceComplexationActivityModel object assuming the presence the Diffuse Double Layer (DDL) model and
/// considering every site as a separate phase.
auto activityModelSurfaceComplexationSiteWithDDL(const SpeciesList& species, ActivityModelSurfaceComplexationSiteParams params) -> ActivityModel
{
    // Initialize the surface
    ComplexationSurface surface = params.surface;

    // The charges of the surface species
    ArrayXd surface_z = surface.charges();

    // Initialize the surface site
    ComplexationSurfaceSite site = surface.sites()[params.site_tag];

    // The indices of site species
    auto indices = site.speciesIndices();

    // The state of the surface and site
    ComplexationSurfaceSiteState site_state;
    ComplexationSurfaceState surface_state;

    // Define the activity model function of the surface complexation site phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the surface complexation
        site_state = site.state(T, P, x);

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        props.ln_a = x.log();

        // If the AqueousPhase has been already evaluated, use the ionic strength to update the electrostatic potential
        if (props.extra["ComplexationSurfaceState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            surface_state = std::any_cast<ComplexationSurfaceState>(props.extra["ComplexationSurfaceState"]);

            // Update surface fractions with calculated site's fractions and surface charge
            surface_state.updateFractions(x, indices);
            surface_state.updateCharge(surface_z);

            // If the AqueousPhase has been already evaluated, use the ionic strength to update the electrostatic potential
            if (props.extra["AqueousMixtureState"].has_value())
            {
                // Export aqueous mixture state via `extra` data member
                const auto& aqueous_state = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

                // Update surface potential using the ionic strength
                surface_state.updatePotential(aqueous_state.Is);
            }
        }
        else
        {
            // Initialize surface state with given temperature and pressure
            surface_state = surface.state(T, P);

            // Update surface fractions with calculated site's fractions and surface charge
            surface_state.updateFractions(x, indices);
            surface_state.updateCharge(surface_z);
        }
        //getchar();

        // Export the surface complexation, site and their states via the `extra` data member
        props.extra["ComplexationSurface"] = surface;
        props.extra["ComplexationSurfaceState"] = surface_state;

        props.extra["ComplexationSurfaceSiteState" + site.name()] = site_state;
        props.extra["ComplexationSurfaceSite" + site.name()] = site;
    };

    return fn;
}

/// Return the SurfaceComplexationActivityModel object assuming the presence the Diffuse Double Layer (DDL) model.
auto activityModelSurfaceComplexationWithDDL(const SpeciesList& species, ActivityModelSurfaceComplexationParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

    // The state of the complexation surface
    ComplexationSurfaceState surface_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the surface complexation
        surface_state = surface.state(T, P, x);

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        props.ln_a = x.log();

        // If the AqueousPhase has been already evaluated, use the ionic strength to update the electrostatic potential
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqueous_state = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Update surface potential using the ionic strength
            surface_state.updatePotential(aqueous_state.Is);
        }

        // Export the surface complexation and its state via the `extra` data member
        props.extra["ComplexationSurfaceState"] = surface_state;
        props.extra["ComplexationSurface"] = surface;
    };

    return fn;
}

/// Return the SurfaceComplexationActivityModel object assuming the presence the Diffuse Double Layer (DDL) model.
auto activityModelSurfaceComplexationWithElectrostatics(const SpeciesList& species, ActivityModelSurfaceComplexationParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

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

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        props.ln_a = x.log();

        // If the AqueousPhase has been already evaluated, use the ionic strength to update the electrostatic potential
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqueous_state = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Update surface potential using the ionic strength
            surface_state.updatePotential(aqueous_state.Is);

            // Update ln of activity and activity coefficients based on the Dzombak & Morel (1990), (2.13), (2.40)
            props.ln_g = z*surface_state.psi*F/(R*T);
            //props.ln_g = 2*z*asinh(surface_state.sigma/(0.1174*sqrt(aqueous_state.Is)));

//            // Update ln of activity coefficients according to the Table 2.12 in Dzombak & Morel (1990)
//            auto eps = 78.4;        // the dielectric constant of water, C/(V*m)
//            auto eps0 = 8.854e-12;  // the permittivity of free space
//            auto kappa = F*sqrt(2*aqueous_state.Is*1e3/(eps*eps0*R*T)); // 1/kappa is the double-layer thickness, m
//            auto psi_GC = surface_state.sigma/(kappa*eps*eps0); // surface potential
//            props.ln_g = z*psi_GC*F/(R*T); // the activity coefficient of a surface species
        }
        props.ln_a += props.ln_g;

        // Export the surface complexation and its state via the `extra` data member
        props.extra["ComplexationSurfaceState"] = surface_state;
        props.extra["ComplexationSurface"] = surface;
    };

    return fn;
}

/// Return the SurfaceComplexationActivityModel object assuming the presence the Diffuse Double Layer (DDL) model and
/// considering every site as a separate phase.
auto activityModelSurfaceComplexationSiteWithElectrostatics(const SpeciesList& species, ActivityModelSurfaceComplexationSiteParams params) -> ActivityModel
{
    // Initialize the surface
    ComplexationSurface surface = params.surface;

    // The charges of the surface species
    ArrayXd surface_z = surface.charges();

    // Initialize the surface site
    ComplexationSurfaceSite site = surface.sites()[params.site_tag];

    // The charges of the site species
    ArrayXd z = site.charges();

    // The indices of site species
    auto indices = site.speciesIndices();

    // The state of the surface and site
    ComplexationSurfaceSiteState site_state;
    ComplexationSurfaceState surface_state;

    // Define the activity model function of the surface complexation site phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the surface complexation
        site_state = site.state(T, P, x);

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        props.ln_a = x.log();

        // If the AqueousPhase has been already evaluated, use the ionic strength to update the electrostatic potential
        if (props.extra["ComplexationSurfaceState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            surface_state = std::any_cast<ComplexationSurfaceState>(props.extra["ComplexationSurfaceState"]);

            // Update surface fractions with calculated site's fractions and surface charge
            surface_state.updateFractions(x, indices);
            surface_state.updateCharge(surface_z);

            // If the AqueousPhase has been already evaluated, use the ionic strength to update the electrostatic potential
            if (props.extra["AqueousMixtureState"].has_value())
            {
                // Export aqueous mixture state via `extra` data member
                const auto& aqueous_state = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

                // Update surface potential using the ionic strength
                surface_state.updatePotential(aqueous_state.Is);

                // Update ln of activity and activity coefficients based on the Dzombak & Morel (1990), (2.13), (2.40)
                props.ln_g = z*surface_state.psi*F/(R*T);
                //props.ln_g = 2*z*asinh(surface_state.sigma/(0.1174*sqrt(aqueous_state.Is)));
                //props.ln_g = 2*asinh(surface_state.sigma/(0.1174*sqrt(aqueous_state.Is)));

//            // Update ln of activity coefficients according to the Table 2.12 in Dzombak & Morel (1990)
//            auto eps = 78.4;        // the dielectric constant of water, C/(V*m)
//            auto eps0 = 8.854e-12;  // the permittivity of free space
//            auto kappa = F*sqrt(2*aqueous_state.Is*1e3/(eps*eps0*R*T)); // 1/kappa is the double-layer thickness, m
//            auto psi_GC = surface_state.sigma/(kappa*eps*eps0); // surface potential
//            props.ln_g = z*psi_GC*F/(R*T); // the activity coefficient of a surface species
            }
            props.ln_a += props.ln_g;

        }
        else
        {
            // Initialize surface state with given temperature and pressure
            surface_state = surface.state(T, P);

            // Update surface fractions with calculated site's fractions and surface charge
                surface_state.updateFractions(x, indices);
                surface_state.updateCharge(surface_z);
        }
        //getchar();

        // Export the surface complexation, site and their states via the `extra` data member
        props.extra["ComplexationSurface"] = surface;
        props.extra["ComplexationSurfaceState"] = surface_state;

        props.extra["ComplexationSurfaceSiteState" + site.name()] = site_state;
        props.extra["ComplexationSurfaceSite" + site.name()] = site;
    };

    return fn;
}

/// Return the Diffuse Double Layer (DDL) activity model based on the Dzombak and Morel (1990) model
auto activityModelDonnanDDL(const SpeciesList& species, ActivityModelDDLParams params) -> ActivityModel
{
    // Create the aqueous ddl_mixture
    AqueousMixture ddl_mixture(species);

    // The array of the ionic species charges
    ArrayXr z = ddl_mixture.charges();

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

        // Add electrostatic correction if the surface complexation activity model has been already considered
        if (props.extra["ComplexationSurfaceState"].has_value())
        {
            // Export surface complexation state via `extra` data member
            auto surface_state = std::any_cast<ComplexationSurfaceState>(props.extra["ComplexationSurfaceState"]);

            // Update activity coefficient using the coulombic correction factor, Langmuir book in (10.34), p. 374
            props.ln_a += -z*F*surface_state.psi/(R*T);
            props.ln_g += -z*F*surface_state.psi/(R*T);
        }
    };

    return fn;
}

/// Return the Diffuse Double Layer (DDL) activity model based on the Dzombak and Morel (1990) model applied on
/// the aqueous solution.
auto activityModelElectrostatics(const SpeciesList& species, ActivityModelDDLParams params) -> ActivityModel
{
    // Create the aqueous ddl_mixture
    AqueousMixture mixture(species);

    // The array of the ionic species charges
    ArrayXr z = mixture.charges();

    // The ddl_state of the aqueous ddl_mixture
    AqueousMixtureState aqueous_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the ddl_state of the aqueous ddl_mixture
        aqueous_state = mixture.state(T, P, x);

        // Export the surface complexation and its ddl_state via the `extra` data member
        props.extra["AqueousMixtureState"] = aqueous_state;

        // Add electrostatic correction if the surface complexation activity model has been already considered
        if (props.extra["ComplexationSurfaceState"].has_value()) {
            // Export surface complexation state via `extra` data member
            auto surface_state = std::any_cast<ComplexationSurfaceState>(props.extra["ComplexationSurfaceState"]);

//            // Update ln of activity coefficients according to the Table 2.12 in Dzombak & Morel (1990)
//            auto eps = 78.5;        // the dielectric constant of water, C/(V*m)
//            auto eps0 = 8.854e-12;  // the permittivity of free space
//            auto kappa = F*sqrt(2*aqueous_state.Is*1e3/(eps*eps0*R*T)); // 1/kappa is the double-layer thickness, m
//            auto psi_GC = surface_state.sigma/(kappa*eps*eps0); // surface potential
//            props.ln_g += -z*psi_GC*F/(R*T); // the activity coefficient of a surface species
//            props.ln_a += -z*psi_GC*F/(R*T);
//
            // Update activity coefficient using the coulombic correction factor, Langmuir book in (10.29) or (10.34)
            props.ln_g += -z*F*surface_state.psi/(R*T);
            props.ln_a += -z*F*surface_state.psi/(R*T);
        }
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

/// Return the SurfaceComplexationActivityModel object assuming no electrostatic effects and considering every site as
/// a separate phase.
auto ActivityModelSurfaceComplexationSiteNoDDL(ActivityModelSurfaceComplexationSiteParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationSiteNoDDL(surface_species, params);
    };
}

/// Return the SurfaceComplexationActivityModel object assuming electrostatic effects.
auto ActivityModelSurfaceComplexationWithDDL(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationWithDDL(surface_species, params);
    };
}

auto ActivityModelSurfaceComplexationWithElectrostatics(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationWithElectrostatics(surface_species, params);
    };
}

auto ActivityModelSurfaceComplexationSiteWithElectrostatics(ActivityModelSurfaceComplexationSiteParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationSiteWithElectrostatics(surface_species, params);
    };
}

/// Return the SurfaceComplexationActivityModel object assuming electrostatic effects and considering every site as
/// a separate phase.
auto ActivityModelSurfaceComplexationSiteWithDDL(ActivityModelSurfaceComplexationSiteParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationSiteWithDDL(surface_species, params);
    };
}

/// Return the Diffuse Double Layer (DDL) activity model based on the Donnan model.
auto ActivityModelDonnanDDL() -> ActivityModelGenerator
{
    ActivityModelDDLParams params;

    return [=](const SpeciesList& species)
    {
        return detail::activityModelDonnanDDL(species, params);
    };
}

/// Return the Diffuse Double Layer (DDL) activity model based on the Donnan model with input params.
auto ActivityModelDonnanDDL(ActivityModelDDLParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return detail::activityModelDonnanDDL(species, params);
    };
}

/// Return the activity model that applied electrostatic effects on aqueous solution.
auto ActivityModelElectrostatics(ActivityModelDDLParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return detail::activityModelElectrostatics(species, params);
    };
}

} // namespace Reaktoro
