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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Thermodynamics/Surface/ActivityModelSurfaceComplexation.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>

using namespace Reaktoro;

/// Initialize mole fractions for the sites' species.
inline auto initializeSurfaceComplexationStrongSiteMoleFractions(const SpeciesList& species) -> ArrayXr
{
    auto idx = [&](auto formula) { return species.indexWithName(formula); };

    ArrayXr n = 1e-6 * ArrayXr::Ones(species.size());
    n[idx("Hfo_sOH")] = 0.1;
    n[idx("Hfo_sOHCa+2" )] = 0.2;
    n[idx("Hfo_sOHSr+2" )] = 0.3;
    return n / n.sum();
}

inline auto initializeSurfaceComplexationWeakSiteMoleFractions(const SpeciesList& species) -> ArrayXr
{
    auto idx = [&](auto formula) { return species.indexWithName(formula); };

    ArrayXr n = 1e-6 * ArrayXr::Ones(species.size());
    n[idx("Hfo_wOH")] = 0.1;
    n[idx("Hfo_wOCa+")] = 0.2;
    n[idx("Hfo_wOSr+")] = 0.3;
    return n / n.sum();
}

/// Return mole fractions for the aqueous species.
inline auto initializeAqueousMoleFractions(const SpeciesList& species) -> ArrayXr
{
    auto idx = [&](auto formula) { return species.indexWithFormula(formula); };

    ArrayXr n(species.size());
    n = 0.1;
    n[idx("H2O")] = 55.508;
    n[idx("H+" )] = 1e-7;
    n[idx("OH-")] = 1e-7;
    n[idx("Na+")] = 0.3;
    n[idx("Cl-")] = 0.3;

    return n / n.sum();
}

namespace test
{
    auto getPhreeqcDatabase(const String& name) -> PhreeqcDatabase;
}

TEST_CASE("Testing ActivityModelSurfaceComplexation", "[ActivityModelSurfaceComplexation]")
{
    // Initialize input data
    const auto T = 25.0;
    const auto P = 1.0;

    // Create custom database
    Database db = test::getPhreeqcDatabase("phreeqc.dat");

    // Create complexation surface
    ComplexationSurface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(60, "m2/g")
        .setMass(4.45, "g");

    // Defined strong site of the complexation surface
    surface_Hfo.addSite("Hfo_w", "_w")
        .setAmount(1e-3, "mol");

    // Defined weak site of the complexation surface
    ComplexationSurfaceSite site_Hfo_s;
    site_Hfo_s.setName("Hfo_s")
        .setAmount(0.025e-3, "mol");
    surface_Hfo.addSite(site_Hfo_s);

    // Fetch all the species with Adsorbed aggregate state
    SpeciesList adsorbed_species = db.species().withAggregateState(AggregateState::Adsorbed);

    // Defined and add surface species
    String selected_absorbed_species = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2 "
                                       "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH";
    String selected_absorbed_species_s = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2";
    String selected_absorbed_species_w = "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH";

    // Select the names of considered absorbed species
    SpeciesList species = adsorbed_species.withNames(selected_absorbed_species);
    surface_Hfo.addSurfaceSpecies(species);

    SpeciesList species_s = adsorbed_species.withNames(selected_absorbed_species_s);
    SpeciesList species_w = adsorbed_species.withNames(selected_absorbed_species_w);

    // The activity model params per site
    ActivityModelSurfaceComplexationSiteParams params_s, params_w;
    params_s.surface = surface_Hfo;
    params_s.site_tag = "_s";
    params_w.surface = surface_Hfo;
    params_w.site_tag = "_w";

    // Initialize corresponding species fractions
    const auto xs = initializeSurfaceComplexationStrongSiteMoleFractions(species_s);
    const auto xw = initializeSurfaceComplexationWeakSiteMoleFractions(species_w);

    SECTION("Checking the activities")
    {
        // --------------------------------------------------------------------
        // Evaluate values for the strong site
        // --------------------------------------------------------------------

        // Construct the activity model function with the given ion exchange species.
        ActivityModel fn = ActivityModelSurfaceComplexationSiteNoDDL(params_s)(species_s);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species_s.size());

        // Evaluate the activity props function
        fn(props, {T, P, xs});

        CHECK( props.ln_a[0]  == Approx(-1.79176)  ); // Hfo_sOH
        CHECK( props.ln_a[1]  == Approx(-1.09862)  ); // Hfo_sOHCa+2
        CHECK( props.ln_a[2]  == Approx(-13.3047)  ); // Hfo_sOH2+
        CHECK( props.ln_a[3]  == Approx(-13.3047)  ); // Hfo_sO-
        CHECK( props.ln_a[4]  == Approx(-0.693151) ); // Hfo_sOHSr+2

        // --------------------------------------------------------------------
        // Evaluate values for the weak site
        // --------------------------------------------------------------------

        // Construct the activity model function with the given ion exchange species.
        fn = ActivityModelSurfaceComplexationSiteNoDDL(params_w)(species_w);

        // Create the ActivityProps object with the results.
        props = ActivityProps::create(species_w.size());

        // Evaluate the activity props function
        fn(props, {T, P, xw});

        CHECK( props.ln_a[0] == Approx(-1.79176)  ); // Hfo_wOH
        CHECK( props.ln_a[1] == Approx(-13.3047)  ); // Hfo_wOH2+
        CHECK( props.ln_a[2] == Approx(-13.3047)  ); // Hfo_wO-
        CHECK( props.ln_a[3] == Approx(-1.09862)  ); // Hfo_wOCa+
        CHECK( props.ln_a[4] == Approx(-0.693152) ); // Hfo_wOSr+
        CHECK( props.ln_a[5] == Approx(-13.3047)  ); // Hfo_wOSrOH
    }

    // Define aqueous species list and corresponding fractions
    const auto species_aq = SpeciesList("H2O H+ OH- Na+ Cl- NaCl");
    const auto x_aq = initializeAqueousMoleFractions(species_aq);

    SECTION("Checking the activities coefficients and activities (calculated based on the parameters fetched from phreeqc.dat)")
    {
        // Create the aqueous mixture
        AqueousMixture mixture(species_aq);

        // The state of the aqueous mixture
        AqueousMixtureState aqstate = mixture.state(T, P, x_aq);

        // --------------------------------------------------------------------
        // Evaluate values for the strong site
        // --------------------------------------------------------------------

        // Construct the activity model function with the given ion exchange species.
        ActivityModel fn = ActivityModelSurfaceComplexationSiteNoDDL(params_s)(species_s);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species_s.size());
        props.extra["AqueousMixtureState"] = aqstate;

        // Evaluate the activity props function
        fn(props, {T, P, xs});

        CHECK( props.ln_a[0]  == Approx(-1.79176)  ); // Hfo_sOH
        CHECK( props.ln_a[1]  == Approx(-1.09862)  ); // Hfo_sOHCa+2
        CHECK( props.ln_a[2]  == Approx(-13.3047)  ); // Hfo_sOH2+
        CHECK( props.ln_a[3]  == Approx(-13.3047)  ); // Hfo_sO-
        CHECK( props.ln_a[4]  == Approx(-0.693151) ); // Hfo_sOHSr+2

        for(auto elem : props.ln_g)
            CHECK( elem  == Approx(0.0) );

        // --------------------------------------------------------------------
        // Evaluate values for the weak site
        // --------------------------------------------------------------------

        // Construct the activity model function with the given ion exchange species.
        fn = ActivityModelSurfaceComplexationSiteNoDDL(params_w)(species_w);

        // Create the ActivityProps object with the results.
        props = ActivityProps::create(species_w.size());
        props.extra["AqueousMixtureState"] = aqstate;

        // Evaluate the activity props function
        fn(props, {T, P, xw});

        CHECK( props.ln_a[0] == Approx(-1.79176)  ); // Hfo_wOH
        CHECK( props.ln_a[1] == Approx(-13.3047)  ); // Hfo_wOH2+
        CHECK( props.ln_a[2] == Approx(-13.3047)  ); // Hfo_wO-
        CHECK( props.ln_a[3] == Approx(-1.09862)  ); // Hfo_wOCa+
        CHECK( props.ln_a[4] == Approx(-0.693152) ); // Hfo_wOSr+
        CHECK( props.ln_a[5] == Approx(-13.3047)  ); // Hfo_wOSrOH
    }

}
