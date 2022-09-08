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
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzerHMW.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>
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

TEST_CASE("Testing ActivityModelPitzerHMW", "[ActivityModelPitzerHMW]")
{
    const auto species = SpeciesList("H2O H+ OH- Na+ Cl- Ca++ HCO3- CO3-- CO2 NaCl HCl NaOH");

    const auto T = 300.0;
    const auto P = 12.3e5;
    const auto x = moleFractions(species);

    // Construct the activity props function with the given aqueous species.
    ActivityModel fn = ActivityModelPitzerHMW()(species);

    // Create the ActivityProps object with the results.
    ActivityProps props = ActivityProps::create(species.size());

    // Evaluate the activity props function
    fn(props, {T, P, x});

    CHECK( exp(props.ln_g[0])  == Approx(1.0024963429) ); // H2O
    CHECK( exp(props.ln_g[1])  == Approx(0.6153629865) ); // H+
    CHECK( exp(props.ln_g[2])  == Approx(0.5225985571) ); // OH-
    CHECK( exp(props.ln_g[3])  == Approx(0.6404881090) ); // Na+
    CHECK( exp(props.ln_g[4])  == Approx(0.6607357243) ); // Cl-
    CHECK( exp(props.ln_g[5])  == Approx(0.1577661672) ); // Ca++
    CHECK( exp(props.ln_g[6])  == Approx(0.6896785630) ); // HCO3-
    CHECK( exp(props.ln_g[7])  == Approx(0.0843383222) ); // CO3--
    CHECK( exp(props.ln_g[8])  == Approx(1.0882824499) ); // CO2
    CHECK( exp(props.ln_g[9])  == Approx(1.0000000000) ); // NaCl
    CHECK( exp(props.ln_g[10]) == Approx(1.0000000000) ); // HCl
    CHECK( exp(props.ln_g[11]) == Approx(1.0000000000) ); // NaOH

    checkActivities(x, props);
}
