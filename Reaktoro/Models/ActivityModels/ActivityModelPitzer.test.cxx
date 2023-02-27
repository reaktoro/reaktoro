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
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzer.hpp>
#include <Reaktoro/Singletons/Elements.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ActivityModelPitzer", "[ActivityModelPitzer]")
{
    PhreeqcDatabase db("pitzer.dat");

    Elements::append(Element("Sg").withMolarMass(0.032066000));

    WHEN("less saline at low temperature - chemical elements are Na, Cl")
    {
        // -----------------------------------------------------------------------------
        // Note: The data below for species names, amounts, and activity coefficients
        // were obtained in this way. The PHREEQC script below was used:
        // ~~~
        // PITZER
        //     -macinnes false
        //
        // SOLUTION
        //     temperature 25.0
        //     units mol/kgw
        //     pH 7.0 charge
        //     Na 0.40
        //     Cl 0.40
        // ~~~
        // Its execution was done using phreeqc script.in script.out
        // path/database/pitzer.dat PHREEQC code was changed to print more
        // precision. The data was collected from script.out using a text editor
        // (Visual Studio Code).
        // -----------------------------------------------------------------------------

        const auto species = SpeciesList("OH- H+ H2O Cl- Na+");

        const auto T = 25.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            1.38935e-07, // OH-
            1.38935e-07, // H+
            5.55062e+01, // H2O
            4.00000e-01, // Cl-
            4.00000e-01, // Na+
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPitzer()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]/ln10  == Approx( -0.177041000 ) ); // OH- (PHREEQC: -0.17703, difference: 0.0062 % )
        CHECK( props.ln_g[1]/ln10  == Approx( -0.109118000 ) ); // H+ (PHREEQC: -0.10911, difference: 0.0073 % )
        CHECK( props.ln_g[2]/ln10  == Approx(  0.000450586 ) ); // H2O (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[3]/ln10  == Approx( -0.159553000 ) ); // Cl- (PHREEQC: -0.15954, difference: 0.0081 % )
        CHECK( props.ln_g[4]/ln10  == Approx( -0.159553000 ) ); // Na+ (PHREEQC: -0.15954, difference: 0.0081 % )
    }

    WHEN("less saline at high temperature - chemical elements are Na, Cl")
    {
        // -----------------------------------------------------------------------------
        // Note: The data below for species names, amounts, and activity coefficients
        // were obtained in this way. The PHREEQC script below was used:
        // ~~~
        // PITZER
        //     -macinnes false
        //
        // SOLUTION
        //     temperature 90.0
        //     units mol/kgw
        //     pH 7.0 charge
        //     Na 0.40
        //     Cl 0.40
        // ~~~
        // Its execution was done using phreeqc script.in script.out
        // path/database/pitzer.dat PHREEQC code was changed to print more
        // precision. The data was collected from script.out using a text editor
        // (Visual Studio Code).
        // -----------------------------------------------------------------------------

        const auto species = SpeciesList("OH- H+ H2O Cl- Na+");

        const auto T = 90.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            9.25551e-07, // OH-
            9.25551e-07, // H+
            5.55062e+01, // H2O
            4.00000e-01, // Cl-
            4.00000e-01, // Na+
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPitzer()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]/ln10  == Approx( -0.195169000 ) ); // OH- (PHREEQC: -0.19486, difference: 0.1586 % )
        CHECK( props.ln_g[1]/ln10  == Approx( -0.148796000 ) ); // H+ (PHREEQC: -0.14848, difference: 0.2128 % )
        CHECK( props.ln_g[2]/ln10  == Approx(  0.000492475 ) ); // H2O (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[3]/ln10  == Approx( -0.178501000 ) ); // Cl- (PHREEQC: -0.17819, difference: 0.1745 % )
        CHECK( props.ln_g[4]/ln10  == Approx( -0.178501000 ) ); // Na+ (PHREEQC: -0.17819, difference: 0.1745 % )
    }

    WHEN("more saline at low temperature - chemical elements are Na, Cl")
    {
        // -----------------------------------------------------------------------------
        // Note: The data below for species names, amounts, and activity coefficients
        // were obtained in this way. The PHREEQC script below was used:
        // ~~~
        // PITZER
        //     -macinnes false
        //
        // SOLUTION
        //     temperature 25.0
        //     units mol/kgw
        //     pH 7.0 charge
        //     Na 4.0
        //     Cl 4.0
        // ~~~
        // Its execution was done using phreeqc script.in script.out
        // path/database/pitzer.dat PHREEQC code was changed to print more
        // precision. The data was collected from script.out using a text editor
        // (Visual Studio Code).
        // -----------------------------------------------------------------------------

        const auto species = SpeciesList("OH- H+ H2O Cl- Na+");

        const auto T = 25.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            8.50571e-08, // OH-
            8.50571e-08, // H+
            5.55062e+01, // H2O
            4.00000e+00, // Cl-
            4.00000e+00, // Na+
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPitzer()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]/ln10  == Approx( -0.272281000 ) ); // OH- (PHREEQC: -0.27226, difference: 0.0077 % )
        CHECK( props.ln_g[1]/ln10  == Approx(  0.348307000 ) ); // H+ (PHREEQC:  0.34830, difference: 0.0020 % )
        CHECK( props.ln_g[2]/ln10  == Approx( -0.011321500 ) ); // H2O (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[3]/ln10  == Approx( -0.106151000 ) ); // Cl- (PHREEQC: -0.10614, difference: 0.0104 % )
        CHECK( props.ln_g[4]/ln10  == Approx( -0.106151000 ) ); // Na+ (PHREEQC: -0.10614, difference: 0.0104 % )
    }

    WHEN("more saline at high temperature - chemical elements are Na, Cl")
    {
        // -----------------------------------------------------------------------------
        // Note: The data below for species names, amounts, and activity coefficients
        // were obtained in this way. The PHREEQC script below was used:
        // ~~~
        // PITZER
        //     -macinnes false
        //
        // SOLUTION
        //     temperature 90.0
        //     units mol/kgw
        //     pH 7.0 charge
        //     Na 4.0
        //     Cl 4.0
        // ~~~
        // Its execution was done using phreeqc script.in script.out
        // path/database/pitzer.dat PHREEQC code was changed to print more
        // precision. The data was collected from script.out using a text editor
        // (Visual Studio Code).
        // -----------------------------------------------------------------------------

        const auto species = SpeciesList("OH- H+ H2O Cl- Na+");

        const auto T = 90.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            6.23937e-07, // OH-
            6.23937e-07, // H+
            5.55062e+01, // H2O
            4.00000e+00, // Cl-
            4.00000e+00, // Na+
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPitzer()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]/ln10  == Approx( -0.286221000 ) ); // OH- (PHREEQC: -0.28559, difference: 0.2209 % )
        CHECK( props.ln_g[1]/ln10  == Approx(  0.220514000 ) ); // H+ (PHREEQC:  0.22112, difference: 0.2741 % )
        CHECK( props.ln_g[2]/ln10  == Approx( -0.010878900 ) ); // H2O (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[3]/ln10  == Approx( -0.125490000 ) ); // Cl- (PHREEQC: -0.12487, difference: 0.4965 % )
        CHECK( props.ln_g[4]/ln10  == Approx( -0.125490000 ) ); // Na+ (PHREEQC: -0.12487, difference: 0.4965 % )
    }

    WHEN("less saline at low temperature - chemical elements are B, Ba, Br, C, Ca, Cl, Fe, K, Li, Mg, Mn, Na, S, Sg, Si, Sr")
    {
        // -----------------------------------------------------------------------------
        // Note: The data below for species names, amounts, and activity coefficients
        // were obtained in this way. The PHREEQC script below was used:
        // ~~~
        // PITZER
        //     -macinnes false
        //
        // SOLUTION
        //     temperature 25.0
        //     units mol/kgw
        //     pH 7.0 charge
        //     B  0.010
        //     Ba 0.120
        //     Br 0.020
        //     C  0.050
        //     Ca 0.050
        //     Cl 0.400
        //     Fe 0.070
        //     K  0.080
        //     Li 0.060
        //     Mg 0.040
        //     Mn 0.060
        //     Na 0.400
        //     S  0.050
        //     Sg 0.010
        //     Si 0.010
        //     Sr 0.140
        // ~~~
        // Its execution was done using phreeqc script.in script.out
        // path/database/pitzer.dat PHREEQC code was changed to print more
        // precision. The data was collected from script.out using a text editor
        // (Visual Studio Code).
        // -----------------------------------------------------------------------------

        const auto species = SpeciesList("OH- H+ H2O B(OH)4- B(OH)3 B3O3(OH)4- B4O5(OH)4-2 Ba+2 Br- CO3-2 HCO3- CO2 Ca+2 CaB(OH)4+ Cl- Fe+2 K+ Li+ MgOH+ Mg+2 MgCO3 MgB(OH)4+ Mn+2 Na+ SO4-2 HSO4- HSg- H2Sg (H2Sg)2 H2SiO4-2 H3SiO4- H4SiO4 Sr+2");

        const auto T = 25.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            8.05577e-01, // OH-
            5.25914e-14, // H+
            5.55062e+01, // H2O
            9.33478e-03, // B(OH)4-
            1.64006e-07, // B(OH)3
            1.41993e-14, // B3O3(OH)4-
            4.06147e-15, // B4O5(OH)4-2
            1.20000e-01, // Ba+2
            2.00000e-02, // Br-
            4.96127e-02, // CO3-2
            1.96310e-06, // HCO3-
            4.98111e-14, // CO2
            4.93890e-02, // Ca+2
            6.10968e-04, // CaB(OH)4+
            4.00000e-01, // Cl-
            7.00000e-02, // Fe+2
            8.00000e-02, // K+
            6.00000e-02, // Li+
            3.46260e-02, // MgOH+
            4.93458e-03, // Mg+2
            3.85342e-04, // MgCO3
            5.40837e-05, // MgB(OH)4+
            6.00000e-02, // Mn+2
            4.00000e-01, // Na+
            5.00000e-02, // SO4-2
            1.02489e-14, // HSO4-
            1.00000e-02, // HSg-
            8.17690e-10, // H2Sg
            4.12410e-20, // (H2Sg)2
            9.79883e-03, // H2SiO4-2
            2.01153e-04, // H3SiO4-
            1.27242e-08, // H4SiO4
            1.40000e-01, // Sr+2
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPitzer()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]/ln10  == Approx( -0.315317000 ) ); // OH- (PHREEQC: -0.31654, difference: 0.3864 % )
        CHECK( props.ln_g[1]/ln10  == Approx( -0.319963000 ) ); // H+ (PHREEQC: -0.31960, difference: 0.1136 % )
        CHECK( props.ln_g[2]/ln10  == Approx(  0.003825070 ) ); // H2O (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[3]/ln10  == Approx( -0.420930000 ) ); // B(OH)4- (PHREEQC: -0.42056, difference: 0.0880 % )
        CHECK( props.ln_g[4]/ln10  == Approx( -0.010631900 ) ); // B(OH)3 (PHREEQC: -0.01063, difference: 0.0179 % )
        CHECK( props.ln_g[5]/ln10  == Approx( -0.440533000 ) ); // B4O5(OH)4-2 (PHREEQC: -0.44017, difference: 0.0825 % )
        CHECK( props.ln_g[6]/ln10  == Approx( -1.686780000 ) ); // B3O3(OH)4- (PHREEQC: -1.68530, difference: 0.0878 % )
        CHECK( props.ln_g[7]/ln10  == Approx( -1.073720000 ) ); // Ba+2 (PHREEQC: -1.07227, difference: 0.1352 % )
        CHECK( props.ln_g[8]/ln10  == Approx( -0.195965000 ) ); // Br- (PHREEQC: -0.19560, difference: 0.1866 % )
        CHECK( props.ln_g[9]/ln10  == Approx( -1.454360000 ) ); // CO3-2 (PHREEQC: -1.45289, difference: 0.1012 % )
        CHECK( props.ln_g[10]/ln10 == Approx( -0.310029000 ) ); // HCO3- (PHREEQC: -0.30966, difference: 0.1192 % )
        CHECK( props.ln_g[11]/ln10 == Approx(  0.043102200 ) ); // CO2 (PHREEQC:  0.04310, difference: 0.0051 % )
        CHECK( props.ln_g[12]/ln10 == Approx( -1.463490000 ) ); // Ca+2 (PHREEQC: -1.49053, difference: 1.8141 % )
        CHECK( props.ln_g[13]/ln10 == Approx( -0.383743000 ) ); // CaB(OH)4+ (PHREEQC: -0.38338, difference: 0.0947 % )
        CHECK( props.ln_g[14]/ln10 == Approx( -0.129128000 ) ); // Cl- (PHREEQC: -0.12908, difference: 0.0372 % )
        CHECK( props.ln_g[15]/ln10 == Approx( -1.337230000 ) ); // Fe+2 (PHREEQC: -1.33576, difference: 0.1100 % )
        CHECK( props.ln_g[16]/ln10 == Approx( -0.236991000 ) ); // K+ (PHREEQC: -0.23663, difference: 0.1526 % )
        CHECK( props.ln_g[17]/ln10 == Approx( -0.302362000 ) ); // Li+ (PHREEQC: -0.30200, difference: 0.1199 % )
        CHECK( props.ln_g[18]/ln10 == Approx( -0.349412000 ) ); // MgOH+ (PHREEQC: -0.34904, difference: 0.1066 % )
        CHECK( props.ln_g[19]/ln10 == Approx( -1.279670000 ) ); // Mg+2 (PHREEQC: -1.27821, difference: 0.1142 % )
        CHECK( props.ln_g[20]/ln10 == Approx(  0.000000000 ) ); // MgCO3 (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[21]/ln10 == Approx( -0.369845000 ) ); // MgB(OH)4+ (PHREEQC: -0.36948, difference: 0.0988 % )
        CHECK( props.ln_g[22]/ln10 == Approx( -1.344300000 ) ); // Mn+2 (PHREEQC: -1.34283, difference: 0.1095 % )
        CHECK( props.ln_g[23]/ln10 == Approx( -0.243301000 ) ); // Na+ (PHREEQC: -0.24294, difference: 0.1486 % )
        CHECK( props.ln_g[24]/ln10 == Approx( -1.323140000 ) ); // SO4-2 (PHREEQC: -1.32168, difference: 0.1105 % )
        CHECK( props.ln_g[25]/ln10 == Approx( -0.253491000 ) ); // HSO4- (PHREEQC: -0.25313, difference: 0.1426 % )
        CHECK( props.ln_g[26]/ln10 == Approx( -0.395627000 ) ); // HSg- (PHREEQC: -0.39526, difference: 0.0929 % )
        CHECK( props.ln_g[27]/ln10 == Approx(  0.035149800 ) ); // H2Sg (PHREEQC:  0.03515, difference: 0.0006 % )
        CHECK( props.ln_g[28]/ln10 == Approx(  0.001912690 ) ); // (H2Sg)2 (PHREEQC:  0.00191, difference: 0.1408 % )
        CHECK( props.ln_g[29]/ln10 == Approx( -1.653040000 ) ); // H2SiO4-2 (PHREEQC: -1.65157, difference: 0.0890 % )
        CHECK( props.ln_g[30]/ln10 == Approx( -0.395627000 ) ); // H3SiO4- (PHREEQC: -0.39526, difference: 0.0929 % )
        CHECK( props.ln_g[31]/ln10 == Approx(  0.036311100 ) ); // H4SiO4 (PHREEQC:  0.03631, difference: 0.0030 % )
        CHECK( props.ln_g[32]/ln10 == Approx( -1.319990000 ) ); // Sr+2 (PHREEQC: -1.31852, difference: 0.1115 % )
    }

    WHEN("less saline at high temperature - chemical elements are B, Ba, Br, C, Ca, Cl, Fe, K, Li, Mg, Mn, Na, S, Sg, Si, Sr")
    {
        // -----------------------------------------------------------------------------
        // Note: The data below for species names, amounts, and activity coefficients
        // were obtained in this way. The PHREEQC script below was used:
        // ~~~
        // PITZER
        //     -macinnes false
        //
        // SOLUTION
        //     temperature 90.0
        //     units mol/kgw
        //     pH 7.0 charge
        //     B  0.010
        //     Ba 0.120
        //     Br 0.020
        //     C  0.050
        //     Ca 0.050
        //     Cl 0.400
        //     Fe 0.070
        //     K  0.080
        //     Li 0.060
        //     Mg 0.040
        //     Mn 0.060
        //     Na 0.400
        //     S  0.050
        //     Sg 0.010
        //     Si 0.010
        //     Sr 0.140
        // ~~~
        // Its execution was done using phreeqc script.in script.out
        // path/database/pitzer.dat PHREEQC code was changed to print more
        // precision. The data was collected from script.out using a text editor
        // (Visual Studio Code).
        // -----------------------------------------------------------------------------

        const auto species = SpeciesList("OH- H+ H2O B(OH)4- B(OH)3 B3O3(OH)4- B4O5(OH)4-2 Ba+2 Br- CO3-2 HCO3- CO2 Ca+2 CaB(OH)4+ Cl- Fe+2 K+ Li+ MgOH+ Mg+2 MgCO3 MgB(OH)4+ Mn+2 Na+ SO4-2 HSO4- HSg- H2Sg (H2Sg)2 H2SiO4-2 H3SiO4- H4SiO4 Sr+2");

        const auto T = 90.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            8.03718e-01, // OH-
            2.66971e-12, // H+
            5.55062e+01, // H2O
            9.58169e-03, // B(OH)4-
            6.36614e-06, // B(OH)3
            2.18584e-11, // B3O3(OH)4-
            8.56360e-12, // B4O5(OH)4-2
            1.20000e-01, // Ba+2
            2.00000e-02, // Br-
            4.97017e-02, // CO3-2
            4.10503e-05, // HCO3-
            4.16122e-11, // CO2
            4.96091e-02, // Ca+2
            3.90859e-04, // CaB(OH)4+
            4.00000e-01, // Cl-
            7.00000e-02, // Fe+2
            8.00000e-02, // K+
            6.00000e-02, // Li+
            3.66429e-02, // MgOH+
            3.07877e-03, // Mg+2
            2.57254e-04, // MgCO3
            2.10797e-05, // MgB(OH)4+
            6.00000e-02, // Mn+2
            4.00000e-01, // Na+
            5.00000e-02, // SO4-2
            2.44360e-12, // HSO4-
            9.99999e-03, // HSg-
            1.17593e-08, // H2Sg
            2.19016e-17, // (H2Sg)2
            9.68619e-03, // H2SiO4-2
            3.13643e-04, // H3SiO4-
            1.64029e-07, // H4SiO4
            1.40000e-01, // Sr+2
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPitzer()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]/ln10  == Approx( -0.361251000 ) ); // OH- (PHREEQC: -0.36188, difference: 0.1738 % )
        CHECK( props.ln_g[1]/ln10  == Approx( -0.389184000 ) ); // H+ (PHREEQC: -0.38817, difference: 0.2612 % )
        CHECK( props.ln_g[2]/ln10  == Approx(  0.004521510 ) ); // H2O (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[3]/ln10  == Approx( -0.480133000 ) ); // B(OH)4- (PHREEQC: -0.47916, difference: 0.2031 % )
        CHECK( props.ln_g[4]/ln10  == Approx( -0.010631900 ) ); // B(OH)3 (PHREEQC: -0.01063, difference: 0.0179 % )
        CHECK( props.ln_g[5]/ln10  == Approx( -0.499812000 ) ); // B4O5(OH)4-2 (PHREEQC: -0.49884, difference: 0.1949 % )
        CHECK( props.ln_g[6]/ln10  == Approx( -1.933250000 ) ); // B3O3(OH)4- (PHREEQC: -1.92926, difference: 0.2068 % )
        CHECK( props.ln_g[7]/ln10  == Approx( -1.297240000 ) ); // Ba+2 (PHREEQC: -1.29335, difference: 0.3008 % )
        CHECK( props.ln_g[8]/ln10  == Approx( -0.215356000 ) ); // Br- (PHREEQC: -0.21439, difference: 0.4506 % )
        CHECK( props.ln_g[9]/ln10  == Approx( -1.641640000 ) ); // CO3-2 (PHREEQC: -1.63765, difference: 0.2436 % )
        CHECK( props.ln_g[10]/ln10 == Approx( -0.371845000 ) ); // HCO3- (PHREEQC: -0.37087, difference: 0.2629 % )
        CHECK( props.ln_g[11]/ln10 == Approx(  0.042842200 ) ); // CO2 (PHREEQC:  0.04284, difference: 0.0051 % )
        CHECK( props.ln_g[12]/ln10 == Approx( -1.670580000 ) ); // Ca+2 (PHREEQC: -1.70216, difference: 1.8553 % )
        CHECK( props.ln_g[13]/ln10 == Approx( -0.447349000 ) ); // CaB(OH)4+ (PHREEQC: -0.44633, difference: 0.2283 % )
        CHECK( props.ln_g[14]/ln10 == Approx( -0.159312000 ) ); // Cl- (PHREEQC: -0.15953, difference: 0.1367 % )
        CHECK( props.ln_g[15]/ln10 == Approx( -1.574740000 ) ); // Fe+2 (PHREEQC: -1.57084, difference: 0.2483 % )
        CHECK( props.ln_g[16]/ln10 == Approx( -0.275637000 ) ); // K+ (PHREEQC: -0.27463, difference: 0.3667 % )
        CHECK( props.ln_g[17]/ln10 == Approx( -0.366432000 ) ); // Li+ (PHREEQC: -0.36542, difference: 0.2769 % )
        CHECK( props.ln_g[18]/ln10 == Approx( -0.412905000 ) ); // MgOH+ (PHREEQC: -0.41189, difference: 0.2464 % )
        CHECK( props.ln_g[19]/ln10 == Approx( -1.502130000 ) ); // Mg+2 (PHREEQC: -1.49823, difference: 0.2603 % )
        CHECK( props.ln_g[20]/ln10 == Approx(  0.000000000 ) ); // MgCO3 (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[21]/ln10 == Approx( -0.433451000 ) ); // MgB(OH)4+ (PHREEQC: -0.43244, difference: 0.2338 % )
        CHECK( props.ln_g[22]/ln10 == Approx( -1.581810000 ) ); // Mn+2 (PHREEQC: -1.57791, difference: 0.2472 % )
        CHECK( props.ln_g[23]/ln10 == Approx( -0.256394000 ) ); // Na+ (PHREEQC: -0.25539, difference: 0.3931 % )
        CHECK( props.ln_g[24]/ln10 == Approx( -1.502600000 ) ); // SO4-2 (PHREEQC: -1.49861, difference: 0.2662 % )
        CHECK( props.ln_g[25]/ln10 == Approx( -0.313764000 ) ); // HSO4- (PHREEQC: -0.31279, difference: 0.3114 % )
        CHECK( props.ln_g[26]/ln10 == Approx( -0.454836000 ) ); // HSg- (PHREEQC: -0.45386, difference: 0.2150 % )
        CHECK( props.ln_g[27]/ln10 == Approx(  0.038135700 ) ); // H2Sg (PHREEQC:  0.03813, difference: 0.0149 % )
        CHECK( props.ln_g[28]/ln10 == Approx(  0.019454800 ) ); // (H2Sg)2 (PHREEQC:  0.01945, difference: 0.0247 % )
        CHECK( props.ln_g[29]/ln10 == Approx( -1.899480000 ) ); // H2SiO4-2 (PHREEQC: -1.89549, difference: 0.2105 % )
        CHECK( props.ln_g[30]/ln10 == Approx( -0.454836000 ) ); // H3SiO4- (PHREEQC: -0.45386, difference: 0.2150 % )
        CHECK( props.ln_g[31]/ln10 == Approx(  0.028118100 ) ); // H4SiO4 (PHREEQC:  0.02812, difference: 0.0068 % )
        CHECK( props.ln_g[32]/ln10 == Approx( -1.521460000 ) ); // Sr+2 (PHREEQC: -1.51756, difference: 0.2570 % )
    }

    WHEN("more saline at low temperature - chemical elements are B, Ba, Br, C, Ca, Cl, Fe, K, Li, Mg, Mn, Na, S, Sg, Si, Sr")
    {
        // -----------------------------------------------------------------------------
        // Note: The data below for species names, amounts, and activity coefficients
        // were obtained in this way. The PHREEQC script below was used:
        // ~~~
        // PITZER
        //     -macinnes false
        //
        // SOLUTION
        //     temperature 25.0
        //     units mol/kgw
        //     pH 7.0 charge
        //     B  0.10
        //     Ba 1.20
        //     Br 0.20
        //     C  0.50
        //     Ca 0.50
        //     Cl 4.00
        //     Fe 0.70
        //     K  0.80
        //     Li 0.60
        //     Mg 0.40
        //     Mn 0.60
        //     Na 4.00
        //     S  0.50
        //     Sg 0.10
        //     Si 0.10
        //     Sr 1.40
        // ~~~
        // Its execution was done using phreeqc script.in script.out
        // path/database/pitzer.dat PHREEQC code was changed to print more
        // precision. The data was collected from script.out using a text editor
        // (Visual Studio Code).
        // -----------------------------------------------------------------------------

        const auto species = SpeciesList("OH- H+ H2O B(OH)4- B(OH)3 B4O5(OH)4-2 B3O3(OH)4- Ba+2 Br- CO3-2 HCO3- CO2 Ca+2 CaB(OH)4+ Cl- Fe+2 K+ Li+ MgOH+ Mg+2 MgCO3 MgB(OH)4+ Mn+2 Na+ SO4-2 HSO4- HSg- H2Sg (H2Sg)2 H2SiO4-2 H3SiO4- H4SiO4 Sr+2");

        const auto T = 25.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            8.00093e+00, // OH-
            1.42500e-15, // H+
            5.55062e+01, // H2O
            9.91916e-02, // B(OH)4-
            7.42402e-08, // B(OH)3
            1.59749e-12, // B4O5(OH)4-2
            6.58421e-14, // B3O3(OH)4-
            1.20000e+00, // Ba+2
            2.00000e-01, // Br-
            4.99860e-01, // CO3-2
            2.60403e-07, // HCO3-
            2.48670e-16, // CO2
            4.99232e-01, // Ca+2
            7.67588e-04, // CaB(OH)4+
            4.00000e+00, // Cl-
            7.00000e-01, // Fe+2
            8.00000e-01, // K+
            6.00000e-01, // Li+
            3.99079e-01, // MgOH+
            7.40542e-04, // Mg+2
            1.40200e-04, // MgCO3
            4.07531e-05, // MgB(OH)4+
            6.00000e-01, // Mn+2
            4.00000e+00, // Na+
            5.00000e-01, // SO4-2
            1.16000e-15, // HSO4-
            1.00000e-01, // HSg-
            1.07989e-10, // H2Sg
            4.16975e-21, // (H2Sg)2
            9.99893e-02, // H2SiO4-2
            1.07469e-05, // H3SiO4-
            1.08392e-11, // H4SiO4
            1.40000e+00, // Sr+2
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPitzer()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]/ln10  == Approx( -0.333687000 ) ); // OH- (PHREEQC: -0.33493, difference: 0.3711 % )
        CHECK( props.ln_g[1]/ln10  == Approx(  0.100480000 ) ); // H+ (PHREEQC:  0.10083, difference: 0.3471 % )
        CHECK( props.ln_g[2]/ln10  == Approx( -0.027975200 ) ); // H2O (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[3]/ln10  == Approx( -0.872613000 ) ); // B(OH)4- (PHREEQC: -0.87223, difference: 0.0439 % )
        CHECK( props.ln_g[4]/ln10  == Approx( -0.070357000 ) ); // B(OH)3 (PHREEQC: -0.07036, difference: 0.0043 % )
        CHECK( props.ln_g[5]/ln10  == Approx( -3.099970000 ) ); // B4O5(OH)4-2 (PHREEQC: -3.09843, difference: 0.0497 % )
        CHECK( props.ln_g[6]/ln10  == Approx( -0.835915000 ) ); // B3O3(OH)4- (PHREEQC: -0.83552, difference: 0.0473 % )
        CHECK( props.ln_g[7]/ln10  == Approx( -1.953580000 ) ); // Ba+2 (PHREEQC: -1.95196, difference: 0.0830 % )
        CHECK( props.ln_g[8]/ln10  == Approx( -0.047441300 ) ); // Br- (PHREEQC: -0.04706, difference: 0.8102 % )
        CHECK( props.ln_g[9]/ln10  == Approx( -2.264760000 ) ); // CO3-2 (PHREEQC: -2.26323, difference: 0.0676 % )
        CHECK( props.ln_g[10]/ln10 == Approx( -0.386501000 ) ); // HCO3- (PHREEQC: -0.38612, difference: 0.0987 % )
        CHECK( props.ln_g[11]/ln10 == Approx(  0.412418000 ) ); // CO2 (PHREEQC:  0.41240, difference: 0.0044 % )
        CHECK( props.ln_g[12]/ln10 == Approx( -3.003410000 ) ); // Ca+2 (PHREEQC: -3.03028, difference: 0.8867 % )
        CHECK( props.ln_g[13]/ln10 == Approx( -0.443235000 ) ); // CaB(OH)4+ (PHREEQC: -0.44286, difference: 0.0847 % )
        CHECK( props.ln_g[14]/ln10 == Approx( -0.179621000 ) ); // Cl- (PHREEQC: -0.17953, difference: 0.0507 % )
        CHECK( props.ln_g[15]/ln10 == Approx( -1.625890000 ) ); // Fe+2 (PHREEQC: -1.62439, difference: 0.0923 % )
        CHECK( props.ln_g[16]/ln10 == Approx(  0.317051000 ) ); // K+ (PHREEQC:  0.31740, difference: 0.1100 % )
        CHECK( props.ln_g[17]/ln10 == Approx(  0.007057200 ) ); // Li+ (PHREEQC:  0.00742, difference: 4.8895 % )
        CHECK( props.ln_g[18]/ln10 == Approx( -1.064470000 ) ); // MgOH+ (PHREEQC: -1.06406, difference: 0.0385 % )
        CHECK( props.ln_g[19]/ln10 == Approx( -1.087990000 ) ); // Mg+2 (PHREEQC: -1.08652, difference: 0.1353 % )
        CHECK( props.ln_g[20]/ln10 == Approx(  0.000000000 ) ); // MgCO3 (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[21]/ln10 == Approx( -0.304255000 ) ); // MgB(OH)4+ (PHREEQC: -0.30388, difference: 0.1234 % )
        CHECK( props.ln_g[22]/ln10 == Approx( -1.898690000 ) ); // Mn+2 (PHREEQC: -1.89717, difference: 0.0801 % )
        CHECK( props.ln_g[23]/ln10 == Approx(  0.429476000 ) ); // Na+ (PHREEQC:  0.42982, difference: 0.0800 % )
        CHECK( props.ln_g[24]/ln10 == Approx( -2.089910000 ) ); // SO4-2 (PHREEQC: -2.08838, difference: 0.0733 % )
        CHECK( props.ln_g[25]/ln10 == Approx( -0.220650000 ) ); // HSO4- (PHREEQC: -0.22028, difference: 0.1680 % )
        CHECK( props.ln_g[26]/ln10 == Approx( -0.766045000 ) ); // HSg- (PHREEQC: -0.76566, difference: 0.0503 % )
        CHECK( props.ln_g[27]/ln10 == Approx(  0.397311000 ) ); // H2Sg (PHREEQC:  0.39729, difference: 0.0053 % )
        CHECK( props.ln_g[28]/ln10 == Approx( -0.037006000 ) ); // (H2Sg)2 (PHREEQC: -0.03700, difference: 0.0162 % )
        CHECK( props.ln_g[29]/ln10 == Approx( -3.157880000 ) ); // H2SiO4-2 (PHREEQC: -3.15632, difference: 0.0494 % )
        CHECK( props.ln_g[30]/ln10 == Approx( -0.766045000 ) ); // H3SiO4- (PHREEQC: -0.76566, difference: 0.0503 % )
        CHECK( props.ln_g[31]/ln10 == Approx(  0.316641000 ) ); // H4SiO4 (PHREEQC:  0.31663, difference: 0.0035 % )
        CHECK( props.ln_g[32]/ln10 == Approx( -1.487300000 ) ); // Sr+2 (PHREEQC: -1.48580, difference: 0.1010 % )
    }

    WHEN("more saline at high temperature - chemical elements are B, Ba, Br, C, Ca, Cl, Fe, K, Li, Mg, Mn, Na, S, Sg, Si, Sr")
    {
        // -----------------------------------------------------------------------------
        // Note: The data below for species names, amounts, and activity coefficients
        // were obtained in this way. The PHREEQC script below was used:
        // ~~~
        // PITZER
        //     -macinnes false
        //
        // SOLUTION
        //     temperature 90.0
        //     units mol/kgw
        //     pH 7.0 charge
        //     B  0.10
        //     Ba 1.20
        //     Br 0.20
        //     C  0.50
        //     Ca 0.50
        //     Cl 4.00
        //     Fe 0.70
        //     K  0.80
        //     Li 0.60
        //     Mg 0.40
        //     Mn 0.60
        //     Na 4.00
        //     S  0.50
        //     Sg 0.10
        //     Si 0.10
        //     Sr 1.40
        // ~~~
        // Its execution was done using phreeqc script.in script.out
        // path/database/pitzer.dat PHREEQC code was changed to print more
        // precision. The data was collected from script.out using a text editor
        // (Visual Studio Code).
        // -----------------------------------------------------------------------------

        const auto species = SpeciesList("OH- H+ H2O B(OH)4- B(OH)3 B4O5(OH)4-2 B3O3(OH)4- Ba+2 Br- CO3-2 HCO3- CO2 Ca+2 CaB(OH)4+ Cl- Fe+2 K+ Li+ MgOH+ Mg+2 MgCO3 MgB(OH)4+ Mn+2 Na+ SO4-2 HSO4- HSg- H2Sg (H2Sg)2 H2SiO4-2 H3SiO4- H4SiO4 Sr+2");

        const auto T = 90.0 + 273.15;
        const auto P = 1.0e+5;

        const auto n = ArrayXr{{
            8.00219e+00, // OH-
            1.86001e-13, // H+
            5.55062e+01, // H2O
            9.97928e-02, // B(OH)4-
            4.13733e-06, // B(OH)3
            4.67224e-09, // B4O5(OH)4-2
            1.41332e-10, // B3O3(OH)4-
            1.20000e+00, // Ba+2
            2.00000e-01, // Br-
            4.99819e-01, // CO3-2
            1.73416e-05, // HCO3-
            1.06493e-12, // CO2
            4.99818e-01, // Ca+2
            1.81582e-04, // CaB(OH)4+
            4.00000e+00, // Cl-
            7.00000e-01, // Fe+2
            8.00000e-01, // K+
            6.00000e-01, // Li+
            3.97855e-01, // MgOH+
            1.95989e-03, // Mg+2
            1.63285e-04, // MgCO3
            2.14796e-05, // MgB(OH)4+
            6.00000e-01, // Mn+2
            4.00000e+00, // Na+
            5.00000e-01, // SO4-2
            1.83087e-13, // HSO4-
            1.00000e-01, // HSg-
            2.41939e-09, // H2Sg
            4.22865e-18, // (H2Sg)2
            9.99784e-02, // H2SiO4-2
            2.16215e-05, // H3SiO4-
            3.50313e-10, // H4SiO4
            1.40000e+00, // Sr+2
        }};

        const auto x = n / n.sum();

        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelPitzer()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( props.ln_g[0]/ln10  == Approx( -0.653708000 ) ); // OH- (PHREEQC: -0.65394, difference: 0.0355 % )
        CHECK( props.ln_g[1]/ln10  == Approx( -0.053114200 ) ); // H+ (PHREEQC: -0.05164, difference: 2.8548 % )
        CHECK( props.ln_g[2]/ln10  == Approx(  0.026530200 ) ); // H2O (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[3]/ln10  == Approx( -1.039030000 ) ); // B(OH)4- (PHREEQC: -1.03764, difference: 0.1340 % )
        CHECK( props.ln_g[4]/ln10  == Approx( -0.070357000 ) ); // B(OH)3 (PHREEQC: -0.07036, difference: 0.0043 % )
        CHECK( props.ln_g[5]/ln10  == Approx( -3.675560000 ) ); // B4O5(OH)4-2 (PHREEQC: -3.66972, difference: 0.1591 % )
        CHECK( props.ln_g[6]/ln10  == Approx( -1.002350000 ) ); // B3O3(OH)4- (PHREEQC: -1.00096, difference: 0.1389 % )
        CHECK( props.ln_g[7]/ln10  == Approx( -1.795770000 ) ); // Ba+2 (PHREEQC: -1.79019, difference: 0.3117 % )
        CHECK( props.ln_g[8]/ln10  == Approx( -0.222314000 ) ); // Br- (PHREEQC: -0.22093, difference: 0.6264 % )
        CHECK( props.ln_g[9]/ln10  == Approx( -2.341480000 ) ); // CO3-2 (PHREEQC: -2.33567, difference: 0.2488 % )
        CHECK( props.ln_g[10]/ln10 == Approx( -0.514021000 ) ); // HCO3- (PHREEQC: -0.51265, difference: 0.2674 % )
        CHECK( props.ln_g[11]/ln10 == Approx(  0.412705000 ) ); // CO2 (PHREEQC:  0.41269, difference: 0.0036 % )
        CHECK( props.ln_g[12]/ln10 == Approx( -3.643610000 ) ); // Ca+2 (PHREEQC: -3.67343, difference: 0.8118 % )
        CHECK( props.ln_g[13]/ln10 == Approx( -0.623739000 ) ); // CaB(OH)4+ (PHREEQC: -0.62223, difference: 0.2425 % )
        CHECK( props.ln_g[14]/ln10 == Approx( -0.130571000 ) ); // Cl- (PHREEQC: -0.13037, difference: 0.1542 % )
        CHECK( props.ln_g[15]/ln10 == Approx( -2.173420000 ) ); // Fe+2 (PHREEQC: -2.16788, difference: 0.2555 % )
        CHECK( props.ln_g[16]/ln10 == Approx(  0.233303000 ) ); // K+ (PHREEQC:  0.23478, difference: 0.6291 % )
        CHECK( props.ln_g[17]/ln10 == Approx( -0.305938000 ) ); // Li+ (PHREEQC: -0.30444, difference: 0.4921 % )
        CHECK( props.ln_g[18]/ln10 == Approx( -1.244940000 ) ); // MgOH+ (PHREEQC: -1.24341, difference: 0.1230 % )
        CHECK( props.ln_g[19]/ln10 == Approx( -1.809450000 ) ); // Mg+2 (PHREEQC: -1.80392, difference: 0.3066 % )
        CHECK( props.ln_g[20]/ln10 == Approx(  0.000000000 ) ); // MgCO3 (PHREEQC:  0.00000, difference: -- % )
        CHECK( props.ln_g[21]/ln10 == Approx( -0.484759000 ) ); // MgB(OH)4+ (PHREEQC: -0.48326, difference: 0.3102 % )
        CHECK( props.ln_g[22]/ln10 == Approx( -2.446250000 ) ); // Mn+2 (PHREEQC: -2.44069, difference: 0.2278 % )
        CHECK( props.ln_g[23]/ln10 == Approx( -0.017460400 ) ); // Na+ (PHREEQC: -0.01595, difference: 9.4696 % )
        CHECK( props.ln_g[24]/ln10 == Approx( -2.881710000 ) ); // SO4-2 (PHREEQC: -2.87584, difference: 0.2041 % )
        CHECK( props.ln_g[25]/ln10 == Approx( -0.386440000 ) ); // HSO4- (PHREEQC: -0.38507, difference: 0.3558 % )
        CHECK( props.ln_g[26]/ln10 == Approx( -0.932493000 ) ); // HSg- (PHREEQC: -0.93111, difference: 0.1485 % )
        CHECK( props.ln_g[27]/ln10 == Approx(  0.427170000 ) ); // H2Sg (PHREEQC:  0.42715, difference: 0.0047 % )
        CHECK( props.ln_g[28]/ln10 == Approx(  0.138415000 ) ); // (H2Sg)2 (PHREEQC:  0.13841, difference: 0.0036 % )
        CHECK( props.ln_g[29]/ln10 == Approx( -3.733470000 ) ); // H2SiO4-2 (PHREEQC: -3.72762, difference: 0.1569 % )
        CHECK( props.ln_g[30]/ln10 == Approx( -0.932493000 ) ); // H3SiO4- (PHREEQC: -0.93111, difference: 0.1485 % )
        CHECK( props.ln_g[31]/ln10 == Approx(  0.239383000 ) ); // H4SiO4 (PHREEQC:  0.23937, difference: 0.0054 % )
        CHECK( props.ln_g[32]/ln10 == Approx( -1.910710000 ) ); // Sr+2 (PHREEQC: -1.90518, difference: 0.2903 % )
    }
}
