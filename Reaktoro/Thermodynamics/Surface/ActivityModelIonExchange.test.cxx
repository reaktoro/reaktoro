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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Thermodynamics/Surface/ActivityModelIonExchange.hpp>
#include <Reaktoro/Thermodynamics/Surface/IonExchangeSurface.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>

using namespace Reaktoro;

/// Initialize mole fractions for the species.
inline auto initializeMoleFractions(const SpeciesList& species) -> ArrayXr
{
    auto idx = [&](auto formula) { return species.indexWithFormula(formula); };

    ArrayXr n = 1e-6 * ArrayXr::Ones(species.size());
    n[idx("X-")] = 0.1;
    n[idx("MgX2")] = 0.1;
    n[idx("CaX2" )] = 0.2;
    n[idx("NaX")] = 0.3;

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

TEST_CASE("Testing ActivityModelIonExchange", "[ActivityModelIonExchange]")
{
    // Load phreeqc database
    PhreeqcDatabase db("phreeqc.dat");

    // Define ion exchange species list
    // Expected species: X- AlOHX2 AlX3 BaX2 CaX2 CdX2 CuX2 FeX2 KX LiX MgX2 MnX2 NH4X NaX PbX2 SrX2 ZnX2
    SpeciesList species = db.species().withAggregateState(AggregateState::IonExchange);

    // Initialize corresponding species fractions
    const auto x = initializeMoleFractions(species);

    SECTION("Checking the species")
    {
        CHECK(species[0].name()  == "X-"     ); // X-
        CHECK(species[1].name()  == "AlOHX2" ); // AlOHX2
        CHECK(species[2].name()  == "AlX3"   ); // AlX3
        CHECK(species[3].name()  == "BaX2"   ); // BaX2
        CHECK(species[4].name()  == "CaX2"   ); // CaX2
        CHECK(species[5].name()  == "CdX2"   ); // CdX2
        CHECK(species[6].name()  == "CuX2"   ); // CuX2
        CHECK(species[7].name()  == "FeX2"   ); // FeX2
        CHECK(species[8].name()  == "KX"     ); // KX
        CHECK(species[9].name()  == "LiX"    ); // LiX
        CHECK(species[10].name() == "MgX2"   ); // MgX2
        CHECK(species[11].name() == "MnX2"   ); // MnX2
        CHECK(species[12].name() == "NH4X"   ); // NH4X
        CHECK(species[13].name() == "NaX"    ); // NaX
        CHECK(species[14].name() == "PbX2"   ); // PbX2
        CHECK(species[15].name() == "SrX2"   ); // SrX2
        CHECK(species[16].name() == "ZnX2"   ); // ZnX2
    }

    SECTION("Checking the charges")
    {
        // Create the aqueous mixture
        IonExchangeSurface surface(species);

        // The numbers of exchanger's equivalents for exchange species
        ArrayXd ze = surface.ze();

        CHECK( ze[0]  == 0 ); // X-
        CHECK( ze[1]  == 2 ); // AlOHX2
        CHECK( ze[2]  == 3 ); // AlX3
        CHECK( ze[3]  == 2 ); // BaX2
        CHECK( ze[4]  == 2 ); // CaX2
        CHECK( ze[5]  == 2 ); // CdX2
        CHECK( ze[6]  == 2 ); // CuX2
        CHECK( ze[7]  == 2 ); // FeX2
        CHECK( ze[8]  == 1 ); // KX
        CHECK( ze[9]  == 1 ); // LiX
        CHECK( ze[10] == 2 ); // MgX2
        CHECK( ze[11] == 2 ); // MnX2
        CHECK( ze[12] == 1 ); // NH4X
        CHECK( ze[13] == 1 ); // NaX
        CHECK( ze[14] == 2 ); // PbX2
        CHECK( ze[15] == 2 ); // SrX2
        CHECK( ze[16] == 2 ); // ZnX2
    }

    // Initialize input data for the ActivityProps
    const auto T = 300.0;
    const auto P = 12.3e5;

    SECTION("Checking the activities")
    {
        // Construct the activity model function with the given ion exchange species.
        ActivityModel fn = ActivityModelIonExchange()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_a[0]  == Approx(  0.0)     ); // X-
        CHECK( props.ln_a[1]  == Approx(-13.017)   ); // AlOHX2
        CHECK( props.ln_a[2]  == Approx(-12.6116)  ); // AlX3
        CHECK( props.ln_a[3]  == Approx(-13.017)   ); // BaX2
        CHECK( props.ln_a[4]  == Approx(-0.810957) ); // CaX2
        CHECK( props.ln_a[5]  == Approx(-13.017)   ); // CdX2
        CHECK( props.ln_a[6]  == Approx(-13.017)   ); // CuX2
        CHECK( props.ln_a[7]  == Approx(-13.017)   ); // FeX2
        CHECK( props.ln_a[8]  == Approx(-13.7102)  ); // KX
        CHECK( props.ln_a[9]  == Approx(-13.7102)  ); // LiX
        CHECK( props.ln_a[10] == Approx(-1.5041)   ); // MgX2
        CHECK( props.ln_a[11] == Approx(-13.017)   ); // MnX2
        CHECK( props.ln_a[12] == Approx(-13.7102)  ); // NH4X
        CHECK( props.ln_a[13] == Approx(-1.09864)  ); // NaX
        CHECK( props.ln_a[14] == Approx(-13.017)   ); // PbX2
        CHECK( props.ln_a[15] == Approx(-13.017)   ); // SrX2
        CHECK( props.ln_a[16] == Approx(-13.017)   ); // ZnX2
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

        // Construct the activity model function with the given ion exchange species.
        ActivityModel fn = ActivityModelIonExchange()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        props.extra["AqueousMixtureState"] = aqstate;

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]  == Approx( 0.00000)  ); // X-
        CHECK( props.ln_g[1]  == Approx(-2.09438)  ); // AlOHX2
        CHECK( props.ln_g[2]  == Approx(-2.31725)  ); // AlX3
        CHECK( props.ln_g[3]  == Approx(-1.47255)  ); // BaX2
        CHECK( props.ln_g[4]  == Approx(-1.29726)  ); // CaX2
        CHECK( props.ln_g[5]  == Approx(-2.09438)  ); // CdX2
        CHECK( props.ln_g[6]  == Approx(-1.31534)  ); // CuX2
        CHECK( props.ln_g[7]  == Approx(-1.31534)  ); // FeX2
        CHECK( props.ln_g[8]  == Approx(-0.41378)  ); // KX
        CHECK( props.ln_g[9]  == Approx(-0.328835) ); // LiX
        CHECK( props.ln_g[10] == Approx(-1.19483)  ); // MgX2
        CHECK( props.ln_g[11] == Approx(-1.31534)  ); // MnX2
        CHECK( props.ln_g[12] == Approx(-0.485979) ); // NH4X
        CHECK( props.ln_g[13] == Approx(-0.324217) ); // NaX
        CHECK( props.ln_g[14] == Approx(-2.09438)  ); // PbX2
        CHECK( props.ln_g[15] == Approx(-1.30042)  ); // SrX2
        CHECK( props.ln_g[16] == Approx(-1.44923)  ); // ZnX2

        CHECK( props.ln_a[0]  == Approx(  0.0000)  ); // X-
        CHECK( props.ln_a[1]  == Approx(-15.1114)  ); // AlOHX2
        CHECK( props.ln_a[2]  == Approx(-14.9288)  ); // AlX3
        CHECK( props.ln_a[3]  == Approx(-14.4896)  ); // BaX2
        CHECK( props.ln_a[4]  == Approx( -2.10821) ); // CaX2
        CHECK( props.ln_a[5]  == Approx(-15.1114)  ); // CdX2
        CHECK( props.ln_a[6]  == Approx(-14.3324)  ); // CuX2
        CHECK( props.ln_a[7]  == Approx(-14.3324)  ); // FeX2
        CHECK( props.ln_a[8]  == Approx(-14.1240)  ); // KX
        CHECK( props.ln_a[9]  == Approx(-14.039)   ); // LiX
        CHECK( props.ln_a[10] == Approx( -2.69894) ); // MgX2
        CHECK( props.ln_a[11] == Approx(-14.3324)  ); // MnX2
        CHECK( props.ln_a[12] == Approx(-14.1962)  ); // NH4X
        CHECK( props.ln_a[13] == Approx(-1.42286)  ); // NaX
        CHECK( props.ln_a[14] == Approx(-15.1114)  ); // PbX2
        CHECK( props.ln_a[15] == Approx(-14.3174)  ); // SrX2
        CHECK( props.ln_a[16] == Approx(-14.4663)  ); // ZnX2
    }
}
