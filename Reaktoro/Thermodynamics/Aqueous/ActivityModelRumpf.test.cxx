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
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelDebyeHuckel.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelRumpf.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
using namespace Reaktoro;

/// Return mole fractions for the species.
inline auto moleFractions(const SpeciesList& species) -> ArrayXr
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

// Check if the activities of the aqueous species are correct assuming activity coefficients are.
inline auto checkActivities(ArrayXrConstRef x, ActivityPropsConstRef props)
{
    const auto iH2O = 0;

    // The concentrations of the species (molalities for solutes, mole fraction for solvent water)
    ArrayXr c = x/(x[iH2O] * waterMolarMass);
    c[iH2O] = x[iH2O];

    for(auto i = 0; i < x.size(); ++i)
    {
        INFO("i = " << i);
        CHECK( exp(props.ln_a[i] - props.ln_g[i]) == Approx(c[i]) );
    }
}

TEST_CASE("Testing ActivityModelRumpf", "[ActivityModelRumpf]")
{
    const auto species = SpeciesList("H2O H+ OH- Na+ Cl- Ca++ HCO3- CO3-- CO2 NaCl HCl NaOH");

    const auto T = 300.0;
    const auto P = 12.3e5;
    const auto x = moleFractions(species);

    // Construct the activity props function with the given aqueous species.
    ActivityModel debyehuckelfn = ActivityModelDebyeHuckel()(species);

    // Create the ActivityProps object with the results.
    ActivityProps props = ActivityProps::create(species.size());

    // Evaluate the activity props function
    debyehuckelfn(props, {T, P, x});

    WHEN("Using ActivityModelRumpf(CO2) constructor")
    {
        ActivityModel fn = ActivityModelRumpf("CO2")(species);
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0])  == Approx(0.927569) ); // H2O
        CHECK( exp(props.ln_g[1])  == Approx(0.725919) ); // H+
        CHECK( exp(props.ln_g[2])  == Approx(0.564033) ); // OH-
        CHECK( exp(props.ln_g[3])  == Approx(0.719478) ); // Na+
        CHECK( exp(props.ln_g[4])  == Approx(0.594087) ); // Cl-
        CHECK( exp(props.ln_g[5])  == Approx(0.244410) ); // Ca++
        CHECK( exp(props.ln_g[6])  == Approx(0.638893) ); // HCO3-
        CHECK( exp(props.ln_g[7])  == Approx(0.178597) ); // CO3--
        CHECK( exp(props.ln_g[8])  == Approx(1.168960) ); // CO2
        CHECK( exp(props.ln_g[9])  == Approx(1.273510) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.273510) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.273510) ); // NaOH

        checkActivities(x, props);
    }

    WHEN("Using ActivityModelRumpf(NaCl) constructor")
    {
        ActivityModel fn = ActivityModelRumpf("NaCl")(species);
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0])  == Approx(0.927569) ); // H2O
        CHECK( exp(props.ln_g[1])  == Approx(0.725919) ); // H+
        CHECK( exp(props.ln_g[2])  == Approx(0.564033) ); // OH-
        CHECK( exp(props.ln_g[3])  == Approx(0.719478) ); // Na+
        CHECK( exp(props.ln_g[4])  == Approx(0.594087) ); // Cl-
        CHECK( exp(props.ln_g[5])  == Approx(0.244410) ); // Ca++
        CHECK( exp(props.ln_g[6])  == Approx(0.638893) ); // HCO3-
        CHECK( exp(props.ln_g[7])  == Approx(0.178597) ); // CO3--
        CHECK( exp(props.ln_g[8])  == Approx(1.273510) ); // CO2
        CHECK( exp(props.ln_g[9])  == Approx(1.168960) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.273510) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.273510) ); // NaOH

        checkActivities(x, props);
    }

    WHEN("A base activity model, such as Debye-Huckel, has not been used previously")
    {
        ActivityModel fn = ActivityModelRumpf("CO2")(species);
        props.extra = {}; // there is no AqueousMixture not AqueousMixtureState in the `extra` data member of `props`

        REQUIRE_THROWS( fn(props, {T, P, x}) );
    }
}
