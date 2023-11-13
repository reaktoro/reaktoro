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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelExtendedUNIQUAC.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ActivityModelExtendedUNIQUAC", "[ActivityModelExtendedUNIQUAC]")
{
    PhreeqcDatabase db("pitzer.dat");

    WHEN("pure water")
    {
        const auto species = SpeciesList("OH- H+ H2O");

        const auto T = 25.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            1.38935e-07, // OH-
            1.38935e-07, // H+
            5.55062e+01, // H2O
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0].val()/ln10  == Approx(-0.0001895461) );  // OH-
        CHECK( props.ln_g[1].val()/ln10  == Approx(-0.0001895777) );  // H+
        CHECK( props.ln_g[2].val()/ln10  == Approx(0.0).scale(1.0) ); // H2O
    }

    // WHEN("pure water - ionic strength is zero")
    // {
    //     const auto species = SpeciesList("OH- H+ H2O");

    //     const auto T = 25.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         1.0e-36,     // OH-
    //         1.0e-36,     // H+
    //         5.55062e+01, // H2O
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx( 0.000000000) ); // OH-
    //     CHECK( props.ln_g[1]/ln10  == Approx( 0.000000000) ); // H+
    //     CHECK( props.ln_g[2]/ln10  == Approx( 0.000000000) ); // H2O
    // }

    // WHEN("less saline at low temperature - chemical elements are Na, Cl")
    // {
    //     // -----------------------------------------------------------------------------
    //     // Note: The data below for species names, amounts, and activity coefficients
    //     // were obtained in this way. The PHREEQC script below was used:
    //     // ~~~
    //     // PITZER
    //     //     -macinnes false
    //     //
    //     // SOLUTION
    //     //     temperature 25.0
    //     //     units mol/kgw
    //     //     pH 7.0 charge
    //     //     Na 0.40
    //     //     Cl 0.40
    //     // ~~~
    //     // Its execution was done using phreeqc script.in script.out
    //     // path/database/pitzer.dat PHREEQC code was changed to print more
    //     // precision. The data was collected from script.out using a text editor
    //     // (Visual Studio Code).
    //     // -----------------------------------------------------------------------------

    //     const auto species = SpeciesList("OH- H+ H2O Cl- Na+");

    //     const auto T = 25.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         1.38935e-07, // OH-
    //         1.38935e-07, // H+
    //         5.55062e+01, // H2O
    //         4.00000e-01, // Cl-
    //         4.00000e-01, // Na+
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx(-0.177033000) ); // OH- (PHREEQC: -0.17703, difference: 1.69e-03 %)
    //     CHECK( props.ln_g[1]/ln10  == Approx(-0.109110000) ); // H+ (PHREEQC: -0.10911, difference: 0.00e+00 %)
    //     CHECK( props.ln_g[2]/ln10  == Approx( 0.000450555) ); // H2O (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[3]/ln10  == Approx(-0.159545000) ); // Cl- (PHREEQC: -0.15954, difference: 3.13e-03 %)
    //     CHECK( props.ln_g[4]/ln10  == Approx(-0.159545000) ); // Na+ (PHREEQC: -0.15954, difference: 3.13e-03 %)
    // }

    // WHEN("less saline at high temperature - chemical elements are Na, Cl")
    // {
    //     // -----------------------------------------------------------------------------
    //     // Note: The data below for species names, amounts, and activity coefficients
    //     // were obtained in this way. The PHREEQC script below was used:
    //     // ~~~
    //     // PITZER
    //     //     -macinnes false
    //     //
    //     // SOLUTION
    //     //     temperature 90.0
    //     //     units mol/kgw
    //     //     pH 7.0 charge
    //     //     Na 0.40
    //     //     Cl 0.40
    //     // ~~~
    //     // Its execution was done using phreeqc script.in script.out
    //     // path/database/pitzer.dat PHREEQC code was changed to print more
    //     // precision. The data was collected from script.out using a text editor
    //     // (Visual Studio Code).
    //     // -----------------------------------------------------------------------------

    //     const auto species = SpeciesList("OH- H+ H2O Cl- Na+");

    //     const auto T = 90.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         9.25551e-07, // OH-
    //         9.25551e-07, // H+
    //         5.55062e+01, // H2O
    //         4.00000e-01, // Cl-
    //         4.00000e-01, // Na+
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx(-0.194858000) ); // OH- (PHREEQC: -0.19486, difference: 1.03e-03 %)
    //     CHECK( props.ln_g[1]/ln10  == Approx(-0.148485000) ); // H+ (PHREEQC: -0.14848, difference: 3.37e-03 %)
    //     CHECK( props.ln_g[2]/ln10  == Approx( 0.000491236) ); // H2O (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[3]/ln10  == Approx(-0.178190000) ); // Cl- (PHREEQC: -0.17819, difference: 0.00e+00 %)
    //     CHECK( props.ln_g[4]/ln10  == Approx(-0.178190000) ); // Na+ (PHREEQC: -0.17819, difference: 0.00e+00 %)
    // }

    // WHEN("more saline at low temperature - chemical elements are Na, Cl")
    // {
    //     // -----------------------------------------------------------------------------
    //     // Note: The data below for species names, amounts, and activity coefficients
    //     // were obtained in this way. The PHREEQC script below was used:
    //     // ~~~
    //     // PITZER
    //     //     -macinnes false
    //     //
    //     // SOLUTION
    //     //     temperature 25.0
    //     //     units mol/kgw
    //     //     pH 7.0 charge
    //     //     Na 4.0
    //     //     Cl 4.0
    //     // ~~~
    //     // Its execution was done using phreeqc script.in script.out
    //     // path/database/pitzer.dat PHREEQC code was changed to print more
    //     // precision. The data was collected from script.out using a text editor
    //     // (Visual Studio Code).
    //     // -----------------------------------------------------------------------------

    //     const auto species = SpeciesList("OH- H+ H2O Cl- Na+");

    //     const auto T = 25.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         8.50571e-08, // OH-
    //         8.50571e-08, // H+
    //         5.55062e+01, // H2O
    //         4.00000e+00, // Cl-
    //         4.00000e+00, // Na+
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx(-0.272265000) ); // OH- (PHREEQC: -0.27226, difference: 1.84e-03 %)
    //     CHECK( props.ln_g[1]/ln10  == Approx( 0.348323000) ); // H+ (PHREEQC:  0.34830, difference: 6.60e-03 %)
    //     CHECK( props.ln_g[2]/ln10  == Approx(-0.011322000) ); // H2O (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[3]/ln10  == Approx(-0.106135000) ); // Cl- (PHREEQC: -0.10614, difference: 4.71e-03 %)
    //     CHECK( props.ln_g[4]/ln10  == Approx(-0.106135000) ); // Na+ (PHREEQC: -0.10614, difference: 4.71e-03 %)
    // }

    // WHEN("more saline at high temperature - chemical elements are Na, Cl")
    // {
    //     // -----------------------------------------------------------------------------
    //     // Note: The data below for species names, amounts, and activity coefficients
    //     // were obtained in this way. The PHREEQC script below was used:
    //     // ~~~
    //     // PITZER
    //     //     -macinnes false
    //     //
    //     // SOLUTION
    //     //     temperature 90.0
    //     //     units mol/kgw
    //     //     pH 7.0 charge
    //     //     Na 4.0
    //     //     Cl 4.0
    //     // ~~~
    //     // Its execution was done using phreeqc script.in script.out
    //     // path/database/pitzer.dat PHREEQC code was changed to print more
    //     // precision. The data was collected from script.out using a text editor
    //     // (Visual Studio Code).
    //     // -----------------------------------------------------------------------------

    //     const auto species = SpeciesList("OH- H+ H2O Cl- Na+");

    //     const auto T = 90.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         6.23937e-07, // OH-
    //         6.23937e-07, // H+
    //         5.55062e+01, // H2O
    //         4.00000e+00, // Cl-
    //         4.00000e+00, // Na+
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx(-0.285593000) ); // OH- (PHREEQC: -0.28559, difference: 1.05e-03 %)
    //     CHECK( props.ln_g[1]/ln10  == Approx( 0.221143000) ); // H+ (PHREEQC:  0.22112, difference: 1.04e-02 %)
    //     CHECK( props.ln_g[2]/ln10  == Approx(-0.010899100) ); // H2O (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[3]/ln10  == Approx(-0.124862000) ); // Cl- (PHREEQC: -0.12487, difference: 6.41e-03 %)
    //     CHECK( props.ln_g[4]/ln10  == Approx(-0.124862000) ); // Na+ (PHREEQC: -0.12487, difference: 6.41e-03 %)
    // }

    // WHEN("less saline at low temperature - chemical elements are B, Ba, Br, C, Ca, Cl, Fe, K, Li, Mg, Mn, Na, S, Sg, Si, Sr")
    // {
    //     // -----------------------------------------------------------------------------
    //     // Note: The data below for species names, amounts, and activity coefficients
    //     // were obtained in this way. The PHREEQC script below was used:
    //     // ~~~
    //     // PITZER
    //     //     -macinnes false
    //     //
    //     // SOLUTION
    //     //     temperature 25.0
    //     //     units mol/kgw
    //     //     pH 7.0 charge
    //     //     B  0.010
    //     //     Ba 0.120
    //     //     Br 0.020
    //     //     C  0.050
    //     //     Ca 0.050
    //     //     Cl 0.400
    //     //     Fe 0.070
    //     //     K  0.080
    //     //     Li 0.060
    //     //     Mg 0.040
    //     //     Mn 0.060
    //     //     Na 0.400
    //     //     S  0.050
    //     //     Sg 0.010
    //     //     Si 0.010
    //     //     Sr 0.140
    //     // ~~~
    //     // Its execution was done using phreeqc script.in script.out
    //     // path/database/pitzer.dat PHREEQC code was changed to print more
    //     // precision. The data was collected from script.out using a text editor
    //     // (Visual Studio Code).
    //     // -----------------------------------------------------------------------------

    //     const auto species = SpeciesList("OH- H+ H2O B(OH)4- B(OH)3 B3O3(OH)4- B4O5(OH)4-2 Ba+2 Br- CO3-2 HCO3- CO2 Ca+2 CaB(OH)4+ Cl- Fe+2 K+ Li+ MgOH+ Mg+2 MgCO3 MgB(OH)4+ Mn+2 Na+ SO4-2 HSO4- HSg- H2Sg (H2Sg)2 H2SiO4-2 H3SiO4- H4SiO4 Sr+2");

    //     const auto T = 25.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         8.05577e-01, // OH-
    //         5.25914e-14, // H+
    //         5.55062e+01, // H2O
    //         9.33478e-03, // B(OH)4-
    //         1.64006e-07, // B(OH)3
    //         1.41993e-14, // B3O3(OH)4-
    //         4.06147e-15, // B4O5(OH)4-2
    //         1.20000e-01, // Ba+2
    //         2.00000e-02, // Br-
    //         4.96127e-02, // CO3-2
    //         1.96310e-06, // HCO3-
    //         4.98111e-14, // CO2
    //         4.93890e-02, // Ca+2
    //         6.10968e-04, // CaB(OH)4+
    //         4.00000e-01, // Cl-
    //         7.00000e-02, // Fe+2
    //         8.00000e-02, // K+
    //         6.00000e-02, // Li+
    //         3.46260e-02, // MgOH+
    //         4.93458e-03, // Mg+2
    //         3.85342e-04, // MgCO3
    //         5.40837e-05, // MgB(OH)4+
    //         6.00000e-02, // Mn+2
    //         4.00000e-01, // Na+
    //         5.00000e-02, // SO4-2
    //         1.02489e-14, // HSO4-
    //         1.00000e-02, // HSg-
    //         8.17690e-10, // H2Sg
    //         4.12410e-20, // (H2Sg)2
    //         9.79883e-03, // H2SiO4-2
    //         2.01153e-04, // H3SiO4-
    //         1.27242e-08, // H4SiO4
    //         1.40000e-01, // Sr+2
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx(-0.315304000) ); // OH- (PHREEQC: -0.31654, difference: 3.90e-01 %)
    //     CHECK( props.ln_g[1]/ln10  == Approx(-0.319948000) ); // H+ (PHREEQC: -0.31960, difference: 1.09e-01 %)
    //     CHECK( props.ln_g[2]/ln10  == Approx( 0.003824830) ); // H2O (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[3]/ln10  == Approx(-0.420916000) ); // B(OH)4- (PHREEQC: -0.42056, difference: 8.46e-02 %)
    //     CHECK( props.ln_g[4]/ln10  == Approx(-0.010631900) ); // B(OH)3 (PHREEQC: -0.01063, difference: 1.79e-02 %)
    //     CHECK( props.ln_g[5]/ln10  == Approx(-0.440519000) ); // B4O5(OH)4-2 (PHREEQC: -0.44017, difference: 7.93e-02 %)
    //     CHECK( props.ln_g[6]/ln10  == Approx(-1.686720000) ); // B3O3(OH)4- (PHREEQC: -1.68530, difference: 8.43e-02 %)
    //     CHECK( props.ln_g[7]/ln10  == Approx(-1.073670000) ); // Ba+2 (PHREEQC: -1.07227, difference: 1.31e-01 %)
    //     CHECK( props.ln_g[8]/ln10  == Approx(-0.195952000) ); // Br- (PHREEQC: -0.19560, difference: 1.80e-01 %)
    //     CHECK( props.ln_g[9]/ln10  == Approx(-1.454300000) ); // CO3-2 (PHREEQC: -1.45289, difference: 9.70e-02 %)
    //     CHECK( props.ln_g[10]/ln10 == Approx(-0.310016000) ); // HCO3- (PHREEQC: -0.30966, difference: 1.15e-01 %)
    //     CHECK( props.ln_g[11]/ln10 == Approx( 0.043102200) ); // CO2 (PHREEQC:  0.04310, difference: 5.10e-03 %)
    //     CHECK( props.ln_g[12]/ln10 == Approx(-1.463430000) ); // Ca+2 (PHREEQC: -1.49053, difference: 1.82e+00 %)
    //     CHECK( props.ln_g[13]/ln10 == Approx(-0.383729000) ); // CaB(OH)4+ (PHREEQC: -0.38338, difference: 9.10e-02 %)
    //     CHECK( props.ln_g[14]/ln10 == Approx(-0.129115000) ); // Cl- (PHREEQC: -0.12908, difference: 2.71e-02 %)
    //     CHECK( props.ln_g[15]/ln10 == Approx(-1.337170000) ); // Fe+2 (PHREEQC: -1.33576, difference: 1.06e-01 %)
    //     CHECK( props.ln_g[16]/ln10 == Approx(-0.236977000) ); // K+ (PHREEQC: -0.23663, difference: 1.47e-01 %)
    //     CHECK( props.ln_g[17]/ln10 == Approx(-0.302347000) ); // Li+ (PHREEQC: -0.30200, difference: 1.15e-01 %)
    //     CHECK( props.ln_g[18]/ln10 == Approx(-0.349397000) ); // MgOH+ (PHREEQC: -0.34904, difference: 1.02e-01 %)
    //     CHECK( props.ln_g[19]/ln10 == Approx(-1.279610000) ); // Mg+2 (PHREEQC: -1.27821, difference: 1.10e-01 %)
    //     CHECK( props.ln_g[20]/ln10 == Approx( 0.000000000) ); // MgCO3 (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[21]/ln10 == Approx(-0.369831000) ); // MgB(OH)4+ (PHREEQC: -0.36948, difference: 9.50e-02 %)
    //     CHECK( props.ln_g[22]/ln10 == Approx(-1.344240000) ); // Mn+2 (PHREEQC: -1.34283, difference: 1.05e-01 %)
    //     CHECK( props.ln_g[23]/ln10 == Approx(-0.243287000) ); // Na+ (PHREEQC: -0.24294, difference: 1.43e-01 %)
    //     CHECK( props.ln_g[24]/ln10 == Approx(-1.323090000) ); // SO4-2 (PHREEQC: -1.32168, difference: 1.07e-01 %)
    //     CHECK( props.ln_g[25]/ln10 == Approx(-0.253477000) ); // HSO4- (PHREEQC: -0.25313, difference: 1.37e-01 %)
    //     CHECK( props.ln_g[26]/ln10 == Approx(-0.395614000) ); // HSg- (PHREEQC: -0.39526, difference: 8.96e-02 %)
    //     CHECK( props.ln_g[27]/ln10 == Approx( 0.035149800) ); // H2Sg (PHREEQC:  0.03515, difference: 5.69e-04 %)
    //     CHECK( props.ln_g[28]/ln10 == Approx( 0.001912690) ); // (H2Sg)2 (PHREEQC:  0.00191, difference: 1.41e-01 %)
    //     CHECK( props.ln_g[29]/ln10 == Approx(-1.652990000) ); // H2SiO4-2 (PHREEQC: -1.65157, difference: 8.60e-02 %)
    //     CHECK( props.ln_g[30]/ln10 == Approx(-0.395614000) ); // H3SiO4- (PHREEQC: -0.39526, difference: 8.96e-02 %)
    //     CHECK( props.ln_g[31]/ln10 == Approx( 0.036311100) ); // H4SiO4 (PHREEQC:  0.03631, difference: 3.03e-03 %)
    //     CHECK( props.ln_g[32]/ln10 == Approx(-1.319930000) ); // Sr+2 (PHREEQC: -1.31852, difference: 1.07e-01 %)
    // }

    // WHEN("less saline at high temperature - chemical elements are B, Ba, Br, C, Ca, Cl, Fe, K, Li, Mg, Mn, Na, S, Sg, Si, Sr")
    // {
    //     // -----------------------------------------------------------------------------
    //     // Note: The data below for species names, amounts, and activity coefficients
    //     // were obtained in this way. The PHREEQC script below was used:
    //     // ~~~
    //     // PITZER
    //     //     -macinnes false
    //     //
    //     // SOLUTION
    //     //     temperature 90.0
    //     //     units mol/kgw
    //     //     pH 7.0 charge
    //     //     B  0.010
    //     //     Ba 0.120
    //     //     Br 0.020
    //     //     C  0.050
    //     //     Ca 0.050
    //     //     Cl 0.400
    //     //     Fe 0.070
    //     //     K  0.080
    //     //     Li 0.060
    //     //     Mg 0.040
    //     //     Mn 0.060
    //     //     Na 0.400
    //     //     S  0.050
    //     //     Sg 0.010
    //     //     Si 0.010
    //     //     Sr 0.140
    //     // ~~~
    //     // Its execution was done using phreeqc script.in script.out
    //     // path/database/pitzer.dat PHREEQC code was changed to print more
    //     // precision. The data was collected from script.out using a text editor
    //     // (Visual Studio Code).
    //     // -----------------------------------------------------------------------------

    //     const auto species = SpeciesList("OH- H+ H2O B(OH)4- B(OH)3 B3O3(OH)4- B4O5(OH)4-2 Ba+2 Br- CO3-2 HCO3- CO2 Ca+2 CaB(OH)4+ Cl- Fe+2 K+ Li+ MgOH+ Mg+2 MgCO3 MgB(OH)4+ Mn+2 Na+ SO4-2 HSO4- HSg- H2Sg (H2Sg)2 H2SiO4-2 H3SiO4- H4SiO4 Sr+2");

    //     const auto T = 90.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         8.03718e-01, // OH-
    //         2.66971e-12, // H+
    //         5.55062e+01, // H2O
    //         9.58169e-03, // B(OH)4-
    //         6.36614e-06, // B(OH)3
    //         2.18584e-11, // B3O3(OH)4-
    //         8.56360e-12, // B4O5(OH)4-2
    //         1.20000e-01, // Ba+2
    //         2.00000e-02, // Br-
    //         4.97017e-02, // CO3-2
    //         4.10503e-05, // HCO3-
    //         4.16122e-11, // CO2
    //         4.96091e-02, // Ca+2
    //         3.90859e-04, // CaB(OH)4+
    //         4.00000e-01, // Cl-
    //         7.00000e-02, // Fe+2
    //         8.00000e-02, // K+
    //         6.00000e-02, // Li+
    //         3.66429e-02, // MgOH+
    //         3.07877e-03, // Mg+2
    //         2.57254e-04, // MgCO3
    //         2.10797e-05, // MgB(OH)4+
    //         6.00000e-02, // Mn+2
    //         4.00000e-01, // Na+
    //         5.00000e-02, // SO4-2
    //         2.44360e-12, // HSO4-
    //         9.99999e-03, // HSg-
    //         1.17593e-08, // H2Sg
    //         2.19016e-17, // (H2Sg)2
    //         9.68619e-03, // H2SiO4-2
    //         3.13643e-04, // H3SiO4-
    //         1.64029e-07, // H4SiO4
    //         1.40000e-01, // Sr+2
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx(-0.360721000) ); // OH- (PHREEQC: -0.36188, difference: 3.20e-01 %)
    //     CHECK( props.ln_g[1]/ln10  == Approx(-0.388612000) ); // H+ (PHREEQC: -0.38817, difference: 1.14e-01 %)
    //     CHECK( props.ln_g[2]/ln10  == Approx( 0.004511950) ); // H2O (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[3]/ln10  == Approx(-0.479603000) ); // B(OH)4- (PHREEQC: -0.47916, difference: 9.25e-02 %)
    //     CHECK( props.ln_g[4]/ln10  == Approx(-0.010631900) ); // B(OH)3 (PHREEQC: -0.01063, difference: 1.79e-02 %)
    //     CHECK( props.ln_g[5]/ln10  == Approx(-0.499283000) ); // B4O5(OH)4-2 (PHREEQC: -0.49884, difference: 8.88e-02 %)
    //     CHECK( props.ln_g[6]/ln10  == Approx(-1.931030000) ); // B3O3(OH)4- (PHREEQC: -1.92926, difference: 9.17e-02 %)
    //     CHECK( props.ln_g[7]/ln10  == Approx(-1.295110000) ); // Ba+2 (PHREEQC: -1.29335, difference: 1.36e-01 %)
    //     CHECK( props.ln_g[8]/ln10  == Approx(-0.214826000) ); // Br- (PHREEQC: -0.21439, difference: 2.03e-01 %)
    //     CHECK( props.ln_g[9]/ln10  == Approx(-1.639410000) ); // CO3-2 (PHREEQC: -1.63765, difference: 1.07e-01 %)
    //     CHECK( props.ln_g[10]/ln10 == Approx(-0.371315000) ); // HCO3- (PHREEQC: -0.37087, difference: 1.20e-01 %)
    //     CHECK( props.ln_g[11]/ln10 == Approx( 0.042842200) ); // CO2 (PHREEQC:  0.04284, difference: 5.14e-03 %)
    //     CHECK( props.ln_g[12]/ln10 == Approx(-1.668450000) ); // Ca+2 (PHREEQC: -1.70216, difference: 1.98e+00 %)
    //     CHECK( props.ln_g[13]/ln10 == Approx(-0.446776000) ); // CaB(OH)4+ (PHREEQC: -0.44633, difference: 9.99e-02 %)
    //     CHECK( props.ln_g[14]/ln10 == Approx(-0.158782000) ); // Cl- (PHREEQC: -0.15953, difference: 4.69e-01 %)
    //     CHECK( props.ln_g[15]/ln10 == Approx(-1.572610000) ); // Fe+2 (PHREEQC: -1.57084, difference: 1.13e-01 %)
    //     CHECK( props.ln_g[16]/ln10 == Approx(-0.275064000) ); // K+ (PHREEQC: -0.27463, difference: 1.58e-01 %)
    //     CHECK( props.ln_g[17]/ln10 == Approx(-0.365859000) ); // Li+ (PHREEQC: -0.36542, difference: 1.20e-01 %)
    //     CHECK( props.ln_g[18]/ln10 == Approx(-0.412333000) ); // MgOH+ (PHREEQC: -0.41189, difference: 1.08e-01 %)
    //     CHECK( props.ln_g[19]/ln10 == Approx(-1.499990000) ); // Mg+2 (PHREEQC: -1.49823, difference: 1.17e-01 %)
    //     CHECK( props.ln_g[20]/ln10 == Approx( 0.000000000) ); // MgCO3 (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[21]/ln10 == Approx(-0.432878000) ); // MgB(OH)4+ (PHREEQC: -0.43244, difference: 1.01e-01 %)
    //     CHECK( props.ln_g[22]/ln10 == Approx(-1.579670000) ); // Mn+2 (PHREEQC: -1.57791, difference: 1.12e-01 %)
    //     CHECK( props.ln_g[23]/ln10 == Approx(-0.255822000) ); // Na+ (PHREEQC: -0.25539, difference: 1.69e-01 %)
    //     CHECK( props.ln_g[24]/ln10 == Approx(-1.500380000) ); // SO4-2 (PHREEQC: -1.49861, difference: 1.18e-01 %)
    //     CHECK( props.ln_g[25]/ln10 == Approx(-0.313235000) ); // HSO4- (PHREEQC: -0.31279, difference: 1.42e-01 %)
    //     CHECK( props.ln_g[26]/ln10 == Approx(-0.454306000) ); // HSg- (PHREEQC: -0.45386, difference: 9.83e-02 %)
    //     CHECK( props.ln_g[27]/ln10 == Approx( 0.038135700) ); // H2Sg (PHREEQC:  0.03813, difference: 1.49e-02 %)
    //     CHECK( props.ln_g[28]/ln10 == Approx( 0.019454800) ); // (H2Sg)2 (PHREEQC:  0.01945, difference: 2.47e-02 %)
    //     CHECK( props.ln_g[29]/ln10 == Approx(-1.897260000) ); // H2SiO4-2 (PHREEQC: -1.89549, difference: 9.34e-02 %)
    //     CHECK( props.ln_g[30]/ln10 == Approx(-0.454306000) ); // H3SiO4- (PHREEQC: -0.45386, difference: 9.83e-02 %)
    //     CHECK( props.ln_g[31]/ln10 == Approx( 0.028118100) ); // H4SiO4 (PHREEQC:  0.02812, difference: 6.76e-03 %)
    //     CHECK( props.ln_g[32]/ln10 == Approx(-1.519320000) ); // Sr+2 (PHREEQC: -1.51756, difference: 1.16e-01 %)
    // }

    // WHEN("more saline at low temperature - chemical elements are B, Ba, Br, C, Ca, Cl, Fe, K, Li, Mg, Mn, Na, S, Sg, Si, Sr")
    // {
    //     // -----------------------------------------------------------------------------
    //     // Note: The data below for species names, amounts, and activity coefficients
    //     // were obtained in this way. The PHREEQC script below was used:
    //     // ~~~
    //     // PITZER
    //     //     -macinnes false
    //     //
    //     // SOLUTION
    //     //     temperature 25.0
    //     //     units mol/kgw
    //     //     pH 7.0 charge
    //     //     B  0.10
    //     //     Ba 1.20
    //     //     Br 0.20
    //     //     C  0.50
    //     //     Ca 0.50
    //     //     Cl 4.00
    //     //     Fe 0.70
    //     //     K  0.80
    //     //     Li 0.60
    //     //     Mg 0.40
    //     //     Mn 0.60
    //     //     Na 4.00
    //     //     S  0.50
    //     //     Sg 0.10
    //     //     Si 0.10
    //     //     Sr 1.40
    //     // ~~~
    //     // Its execution was done using phreeqc script.in script.out
    //     // path/database/pitzer.dat PHREEQC code was changed to print more
    //     // precision. The data was collected from script.out using a text editor
    //     // (Visual Studio Code).
    //     // -----------------------------------------------------------------------------

    //     const auto species = SpeciesList("OH- H+ H2O B(OH)4- B(OH)3 B4O5(OH)4-2 B3O3(OH)4- Ba+2 Br- CO3-2 HCO3- CO2 Ca+2 CaB(OH)4+ Cl- Fe+2 K+ Li+ MgOH+ Mg+2 MgCO3 MgB(OH)4+ Mn+2 Na+ SO4-2 HSO4- HSg- H2Sg (H2Sg)2 H2SiO4-2 H3SiO4- H4SiO4 Sr+2");

    //     const auto T = 25.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         8.00093e+00, // OH-
    //         1.42500e-15, // H+
    //         5.55062e+01, // H2O
    //         9.91916e-02, // B(OH)4-
    //         7.42402e-08, // B(OH)3
    //         1.59749e-12, // B4O5(OH)4-2
    //         6.58421e-14, // B3O3(OH)4-
    //         1.20000e+00, // Ba+2
    //         2.00000e-01, // Br-
    //         4.99860e-01, // CO3-2
    //         2.60403e-07, // HCO3-
    //         2.48670e-16, // CO2
    //         4.99232e-01, // Ca+2
    //         7.67588e-04, // CaB(OH)4+
    //         4.00000e+00, // Cl-
    //         7.00000e-01, // Fe+2
    //         8.00000e-01, // K+
    //         6.00000e-01, // Li+
    //         3.99079e-01, // MgOH+
    //         7.40542e-04, // Mg+2
    //         1.40200e-04, // MgCO3
    //         4.07531e-05, // MgB(OH)4+
    //         6.00000e-01, // Mn+2
    //         4.00000e+00, // Na+
    //         5.00000e-01, // SO4-2
    //         1.16000e-15, // HSO4-
    //         1.00000e-01, // HSg-
    //         1.07989e-10, // H2Sg
    //         4.16975e-21, // (H2Sg)2
    //         9.99893e-02, // H2SiO4-2
    //         1.07469e-05, // H3SiO4-
    //         1.08392e-11, // H4SiO4
    //         1.40000e+00, // Sr+2
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx(-0.333664000) ); // OH- (PHREEQC: -0.33493, difference: 3.78e-01 %)
    //     CHECK( props.ln_g[1]/ln10  == Approx( 0.100507000) ); // H+ (PHREEQC:  0.10083, difference: 3.20e-01 %)
    //     CHECK( props.ln_g[2]/ln10  == Approx(-0.027978700) ); // H2O (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[3]/ln10  == Approx(-0.872589000) ); // B(OH)4- (PHREEQC: -0.87223, difference: 4.12e-02 %)
    //     CHECK( props.ln_g[4]/ln10  == Approx(-0.070357000) ); // B(OH)3 (PHREEQC: -0.07036, difference: 4.26e-03 %)
    //     CHECK( props.ln_g[5]/ln10  == Approx(-3.099870000) ); // B4O5(OH)4-2 (PHREEQC: -3.09843, difference: 4.65e-02 %)
    //     CHECK( props.ln_g[6]/ln10  == Approx(-0.835892000) ); // B3O3(OH)4- (PHREEQC: -0.83552, difference: 4.45e-02 %)
    //     CHECK( props.ln_g[7]/ln10  == Approx(-1.953490000) ); // Ba+2 (PHREEQC: -1.95196, difference: 7.84e-02 %)
    //     CHECK( props.ln_g[8]/ln10  == Approx(-0.047417600) ); // Br- (PHREEQC: -0.04706, difference: 7.60e-01 %)
    //     CHECK( props.ln_g[9]/ln10  == Approx(-2.264660000) ); // CO3-2 (PHREEQC: -2.26323, difference: 6.32e-02 %)
    //     CHECK( props.ln_g[10]/ln10 == Approx(-0.386477000) ); // HCO3- (PHREEQC: -0.38612, difference: 9.25e-02 %)
    //     CHECK( props.ln_g[11]/ln10 == Approx( 0.412418000) ); // CO2 (PHREEQC:  0.41240, difference: 4.36e-03 %)
    //     CHECK( props.ln_g[12]/ln10 == Approx(-3.003310000) ); // Ca+2 (PHREEQC: -3.03028, difference: 8.90e-01 %)
    //     CHECK( props.ln_g[13]/ln10 == Approx(-0.443207000) ); // CaB(OH)4+ (PHREEQC: -0.44286, difference: 7.84e-02 %)
    //     CHECK( props.ln_g[14]/ln10 == Approx(-0.179597000) ); // Cl- (PHREEQC: -0.17953, difference: 3.73e-02 %)
    //     CHECK( props.ln_g[15]/ln10 == Approx(-1.625790000) ); // Fe+2 (PHREEQC: -1.62439, difference: 8.62e-02 %)
    //     CHECK( props.ln_g[16]/ln10 == Approx( 0.317079000) ); // K+ (PHREEQC:  0.31740, difference: 1.01e-01 %)
    //     CHECK( props.ln_g[17]/ln10 == Approx( 0.007083930) ); // Li+ (PHREEQC:  0.00742, difference: 4.52e+00 %)
    //     CHECK( props.ln_g[18]/ln10 == Approx(-1.064450000) ); // MgOH+ (PHREEQC: -1.06406, difference: 3.67e-02 %)
    //     CHECK( props.ln_g[19]/ln10 == Approx(-1.087900000) ); // Mg+2 (PHREEQC: -1.08652, difference: 1.27e-01 %)
    //     CHECK( props.ln_g[20]/ln10 == Approx( 0.000000000) ); // MgCO3 (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[21]/ln10 == Approx(-0.304228000) ); // MgB(OH)4+ (PHREEQC: -0.30388, difference: 1.15e-01 %)
    //     CHECK( props.ln_g[22]/ln10 == Approx(-1.898600000) ); // Mn+2 (PHREEQC: -1.89717, difference: 7.54e-02 %)
    //     CHECK( props.ln_g[23]/ln10 == Approx( 0.429503000) ); // Na+ (PHREEQC:  0.42982, difference: 7.38e-02 %)
    //     CHECK( props.ln_g[24]/ln10 == Approx(-2.089810000) ); // SO4-2 (PHREEQC: -2.08838, difference: 6.85e-02 %)
    //     CHECK( props.ln_g[25]/ln10 == Approx(-0.220626000) ); // HSO4- (PHREEQC: -0.22028, difference: 1.57e-01 %)
    //     CHECK( props.ln_g[26]/ln10 == Approx(-0.766021000) ); // HSg- (PHREEQC: -0.76566, difference: 4.71e-02 %)
    //     CHECK( props.ln_g[27]/ln10 == Approx( 0.397311000) ); // H2Sg (PHREEQC:  0.39729, difference: 5.29e-03 %)
    //     CHECK( props.ln_g[28]/ln10 == Approx(-0.037006000) ); // (H2Sg)2 (PHREEQC: -0.03700, difference: 1.62e-02 %)
    //     CHECK( props.ln_g[29]/ln10 == Approx(-3.157770000) ); // H2SiO4-2 (PHREEQC: -3.15632, difference: 4.59e-02 %)
    //     CHECK( props.ln_g[30]/ln10 == Approx(-0.766021000) ); // H3SiO4- (PHREEQC: -0.76566, difference: 4.71e-02 %)
    //     CHECK( props.ln_g[31]/ln10 == Approx( 0.316641000) ); // H4SiO4 (PHREEQC:  0.31663, difference: 3.47e-03 %)
    //     CHECK( props.ln_g[32]/ln10 == Approx(-1.487200000) ); // Sr+2 (PHREEQC: -1.48580, difference: 9.42e-02 %)
    // }

    // WHEN("more saline at high temperature - chemical elements are B, Ba, Br, C, Ca, Cl, Fe, K, Li, Mg, Mn, Na, S, Sg, Si, Sr")
    // {
    //     // -----------------------------------------------------------------------------
    //     // Note: The data below for species names, amounts, and activity coefficients
    //     // were obtained in this way. The PHREEQC script below was used:
    //     // ~~~
    //     // PITZER
    //     //     -macinnes false
    //     //
    //     // SOLUTION
    //     //     temperature 90.0
    //     //     units mol/kgw
    //     //     pH 7.0 charge
    //     //     B  0.10
    //     //     Ba 1.20
    //     //     Br 0.20
    //     //     C  0.50
    //     //     Ca 0.50
    //     //     Cl 4.00
    //     //     Fe 0.70
    //     //     K  0.80
    //     //     Li 0.60
    //     //     Mg 0.40
    //     //     Mn 0.60
    //     //     Na 4.00
    //     //     S  0.50
    //     //     Sg 0.10
    //     //     Si 0.10
    //     //     Sr 1.40
    //     // ~~~
    //     // Its execution was done using phreeqc script.in script.out
    //     // path/database/pitzer.dat PHREEQC code was changed to print more
    //     // precision. The data was collected from script.out using a text editor
    //     // (Visual Studio Code).
    //     // -----------------------------------------------------------------------------

    //     const auto species = SpeciesList("OH- H+ H2O B(OH)4- B(OH)3 B4O5(OH)4-2 B3O3(OH)4- Ba+2 Br- CO3-2 HCO3- CO2 Ca+2 CaB(OH)4+ Cl- Fe+2 K+ Li+ MgOH+ Mg+2 MgCO3 MgB(OH)4+ Mn+2 Na+ SO4-2 HSO4- HSg- H2Sg (H2Sg)2 H2SiO4-2 H3SiO4- H4SiO4 Sr+2");

    //     const auto T = 90.0 + 273.15;
    //     const auto P = 1.0e+5;

    //     const auto n = ArrayXr{{
    //         8.00219e+00, // OH-
    //         1.86001e-13, // H+
    //         5.55062e+01, // H2O
    //         9.97928e-02, // B(OH)4-
    //         4.13733e-06, // B(OH)3
    //         4.67224e-09, // B4O5(OH)4-2
    //         1.41332e-10, // B3O3(OH)4-
    //         1.20000e+00, // Ba+2
    //         2.00000e-01, // Br-
    //         4.99819e-01, // CO3-2
    //         1.73416e-05, // HCO3-
    //         1.06493e-12, // CO2
    //         4.99818e-01, // Ca+2
    //         1.81582e-04, // CaB(OH)4+
    //         4.00000e+00, // Cl-
    //         7.00000e-01, // Fe+2
    //         8.00000e-01, // K+
    //         6.00000e-01, // Li+
    //         3.97855e-01, // MgOH+
    //         1.95989e-03, // Mg+2
    //         1.63285e-04, // MgCO3
    //         2.14796e-05, // MgB(OH)4+
    //         6.00000e-01, // Mn+2
    //         4.00000e+00, // Na+
    //         5.00000e-01, // SO4-2
    //         1.83087e-13, // HSO4-
    //         1.00000e-01, // HSg-
    //         2.41939e-09, // H2Sg
    //         4.22865e-18, // (H2Sg)2
    //         9.99784e-02, // H2SiO4-2
    //         2.16215e-05, // H3SiO4-
    //         3.50313e-10, // H4SiO4
    //         1.40000e+00, // Sr+2
    //     }};

    //     const auto x = n / n.sum();

    //     // Construct the activity props function with the given aqueous species.
    //     ActivityModel fn = ActivityModelExtendedUNIQUAC()(species);

    //     // Create the ActivityProps object with the results.
    //     ActivityProps props = ActivityProps::create(species.size());

    //     // Evaluate the activity props function
    //     fn(props, {T, P, x});

    //     CHECK( props.ln_g[0]/ln10  == Approx(-0.652777000) ); // OH- (PHREEQC: -0.65394, difference: 1.78e-01 %)
    //     CHECK( props.ln_g[1]/ln10  == Approx(-0.052050400) ); // H+ (PHREEQC: -0.05164, difference: 7.95e-01 %)
    //     CHECK( props.ln_g[2]/ln10  == Approx( 0.026394400) ); // H2O (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[3]/ln10  == Approx(-1.038100000) ); // B(OH)4- (PHREEQC: -1.03764, difference: 4.43e-02 %)
    //     CHECK( props.ln_g[4]/ln10  == Approx(-0.070357000) ); // B(OH)3 (PHREEQC: -0.07036, difference: 4.26e-03 %)
    //     CHECK( props.ln_g[5]/ln10  == Approx(-3.671520000) ); // B4O5(OH)4-2 (PHREEQC: -3.66972, difference: 4.91e-02 %)
    //     CHECK( props.ln_g[6]/ln10  == Approx(-1.001420000) ); // B3O3(OH)4- (PHREEQC: -1.00096, difference: 4.60e-02 %)
    //     CHECK( props.ln_g[7]/ln10  == Approx(-1.791990000) ); // Ba+2 (PHREEQC: -1.79019, difference: 1.01e-01 %)
    //     CHECK( props.ln_g[8]/ln10  == Approx(-0.221383000) ); // Br- (PHREEQC: -0.22093, difference: 2.05e-01 %)
    //     CHECK( props.ln_g[9]/ln10  == Approx(-2.337440000) ); // CO3-2 (PHREEQC: -2.33567, difference: 7.58e-02 %)
    //     CHECK( props.ln_g[10]/ln10 == Approx(-0.513090000) ); // HCO3- (PHREEQC: -0.51265, difference: 8.58e-02 %)
    //     CHECK( props.ln_g[11]/ln10 == Approx( 0.412705000) ); // CO2 (PHREEQC:  0.41269, difference: 3.63e-03 %)
    //     CHECK( props.ln_g[12]/ln10 == Approx(-3.639830000) ); // Ca+2 (PHREEQC: -3.67343, difference: 9.15e-01 %)
    //     CHECK( props.ln_g[13]/ln10 == Approx(-0.622675000) ); // CaB(OH)4+ (PHREEQC: -0.62223, difference: 7.15e-02 %)
    //     CHECK( props.ln_g[14]/ln10 == Approx(-0.129640000) ); // Cl- (PHREEQC: -0.13037, difference: 5.60e-01 %)
    //     CHECK( props.ln_g[15]/ln10 == Approx(-2.169640000) ); // Fe+2 (PHREEQC: -2.16788, difference: 8.12e-02 %)
    //     CHECK( props.ln_g[16]/ln10 == Approx( 0.234367000) ); // K+ (PHREEQC:  0.23478, difference: 1.76e-01 %)
    //     CHECK( props.ln_g[17]/ln10 == Approx(-0.304874000) ); // Li+ (PHREEQC: -0.30444, difference: 1.43e-01 %)
    //     CHECK( props.ln_g[18]/ln10 == Approx(-1.243880000) ); // MgOH+ (PHREEQC: -1.24341, difference: 3.78e-02 %)
    //     CHECK( props.ln_g[19]/ln10 == Approx(-1.805670000) ); // Mg+2 (PHREEQC: -1.80392, difference: 9.70e-02 %)
    //     CHECK( props.ln_g[20]/ln10 == Approx( 0.000000000) ); // MgCO3 (PHREEQC:  0.00000, difference: -- %)
    //     CHECK( props.ln_g[21]/ln10 == Approx(-0.483695000) ); // MgB(OH)4+ (PHREEQC: -0.48326, difference: 9.00e-02 %)
    //     CHECK( props.ln_g[22]/ln10 == Approx(-2.442470000) ); // Mn+2 (PHREEQC: -2.44069, difference: 7.29e-02 %)
    //     CHECK( props.ln_g[23]/ln10 == Approx(-0.016396600) ); // Na+ (PHREEQC: -0.01595, difference: 2.80e+00 %)
    //     CHECK( props.ln_g[24]/ln10 == Approx(-2.877660000) ); // SO4-2 (PHREEQC: -2.87584, difference: 6.33e-02 %)
    //     CHECK( props.ln_g[25]/ln10 == Approx(-0.385508000) ); // HSO4- (PHREEQC: -0.38507, difference: 1.14e-01 %)
    //     CHECK( props.ln_g[26]/ln10 == Approx(-0.931561000) ); // HSg- (PHREEQC: -0.93111, difference: 4.84e-02 %)
    //     CHECK( props.ln_g[27]/ln10 == Approx( 0.427170000) ); // H2Sg (PHREEQC:  0.42715, difference: 4.68e-03 %)
    //     CHECK( props.ln_g[28]/ln10 == Approx( 0.138415000) ); // (H2Sg)2 (PHREEQC:  0.13841, difference: 3.61e-03 %)
    //     CHECK( props.ln_g[29]/ln10 == Approx(-3.729420000) ); // H2SiO4-2 (PHREEQC: -3.72762, difference: 4.83e-02 %)
    //     CHECK( props.ln_g[30]/ln10 == Approx(-0.931561000) ); // H3SiO4- (PHREEQC: -0.93111, difference: 4.84e-02 %)
    //     CHECK( props.ln_g[31]/ln10 == Approx( 0.239383000) ); // H4SiO4 (PHREEQC:  0.23937, difference: 5.43e-03 %)
    //     CHECK( props.ln_g[32]/ln10 == Approx(-1.906930000) ); // Sr+2 (PHREEQC: -1.90518, difference: 9.19e-02 %)
    // }
}
