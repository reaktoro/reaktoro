// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2021 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// // Catch includes
// #include <catch2/catch.hpp>
// using namespace Catch;

// // C++ includes
// #include <limits>

// // Supkrit includes
// #include <Supkrit/Params.hpp>
// #include <Supkrit/SerializationYAML.hpp>
// using namespace Atomik;
// using namespace Supkrit;

// std::string params_mk = R"(
// Gf: -3679250.6
// Hf: -3876463.4
// Sr: 209.32552
// Vr: 9.281e-05
// a: 251.41656
// b: 0.0476976
// c: -4769760.0
// Tmax: 1700.0
// )";

// std::string params_hkfmk1 = R"(
// Gf: -3708312.7
// Hf: -3931621.1
// Sr: 207.14984
// Vr: 0.00010025
// a: [258.1528, 342.58592]
// b: [0.0581576, 0.014869936]
// c: [-6280184.0, -20984434.0]
// Ttr: [473.0]
// Htr: [.nan]
// Vtr: [.nan]
// dPdTtr: [.nan]
// Tmax: 1200.0
// )";

// std::string params_hkfmk2 = R"(
// Gf: -39522.064
// Hf: -31589.2
// Sr: 143.5112
// Vr: 3.42e-05
// a: [65.39592, 7.610696, 90.3744]
// b: [0.0359824, 0.221752, 0.0]
// c: [0.0, 0.0, 0.0]
// Ttr: [450.0, 620.0]
// Htr: [3974.8, 2510.4]
// Vtr: [.nan, .nan]
// dPdTtr: [.nan, .nan]
// Tmax: 1000.0
// )";

// std::string params_hkfmk3 = R"(
// Gf: .nan
// Hf: .nan
// Sr: 286.604
// Vr: 0.0001432
// a: [369.61456, 409.65544, 488.18912, 461.24416]
// b: [0.22643808, 0.14786256, 0.03112896, 0.04217472]
// c: [-7556304.0, -4167264.0, -1937192.0, -1937192.0]
// Ttr: [848.0, 950.0, 1050.0]
// Htr: [.nan, .nan, .nan]
// Vtr: [.nan, .nan, .nan]
// dPdTtr: [.nan, .nan, .nan]
// Tmax: 1100.0
// )";

// std::string params_whkf = R"(
// Ttr: 273.16
// Str: 63.312288
// Gtr: -235517.36
// Htr: -287721.13
// Utr: -284039.21
// Atr: -231856.36
// )";

// std::string params_hkf = R"(
// Gf: 39371.44
// Hf: -151084.24
// Sr: 197.4848
// a1: 5.8268894e-05
// a2: 8251.0572
// a3: 0.00049988758
// a4: -150377.14
// c1: 384.55521
// c2: 116047.42
// wref: -156816.32
// charge: 0.0
// )";

// std::string params_hp = R"(
// Gf: -4937500.0
// Hf: -5260650.0
// Sr: 342.0
// Vr: 0.00011525
// a: 677.3
// b: 0.0
// c: -3772700.0
// d: -5044.0
// alpha: 2.12e-05
// kappa: 190000000000.0
// kappap: 2.98
// kappapp: -1.6e-11
// numatoms: 20.0
// )";

// std::string params_hpg = R"(
// Gf: -50710.0
// Hf: -74810.0
// Sr: 186.26
// Vr: 0.0
// a: 150.1
// b: 0.002063
// c: 3427700.0
// d: -2650.4
// )";

// std::string params_hpl = R"(
// Gf: -2192340.0
// Hf: -2307040.0
// Sr: 127.6
// Vr: 5.16e-05
// a: 247.5
// b: -0.003206
// c: 0.0
// d: -2051.9
// alpha: 2.9e-05
// kappa: 98500000000.0
// kappap: 4.07
// kappapp: -4.1e-11
// numatoms: 7.0
// Tcr: 1710.0
// Smax: 10.03
// Vmax: 5.0e-07
// )";

// using vec = std::vector<double>;

// inline auto allNaN(const vec& v)
// {
//     for(auto&& x : v)
//         if(!std::isnan(x)) return false;
//     return true;
// }

// TEST_CASE("Testing Serialization", "[Serialization]")
// {
//     const double nan = std::numeric_limits<double>::quiet_NaN();

//     SECTION("Testing YAML serialization of ParamsMaierKelley")
//     {
//         ParamsMaierKelley params = yaml(params_mk);
//         CHECK( params.Gf == -3679250.6 );
//         CHECK( params.Hf == -3876463.4 );
//         CHECK( params.Sr == 209.32552 );
//         CHECK( params.Vr == 9.281e-05 );
//         CHECK( params.a == 251.41656 );
//         CHECK( params.b == 0.0476976 );
//         CHECK( params.c == -4769760.0 );
//         CHECK( params.Tmax == 1700.0 );
//     }

//     SECTION("Testing YAML serialization of ParamsMaierKelleyHKF")
//     {
//         ParamsMaierKelleyHKF params = yaml(params_hkfmk1);
//         CHECK( params.Gf == -3708312.7 );
//         CHECK( params.Hf == -3931621.1 );
//         CHECK( params.Sr == 207.14984 );
//         CHECK( params.Vr == 0.00010025 );
//         CHECK( params.a == vec{ 258.1528, 342.58592 } );
//         CHECK( params.b == vec{ 0.0581576, 0.014869936 } );
//         CHECK( params.c == vec{ -6280184.0, -20984434.0 } );
//         CHECK( params.Ttr == vec{ 473.0 } );
//         CHECK( allNaN( params.Htr ) );
//         CHECK( allNaN( params.Vtr ) );
//         CHECK( allNaN( params.dPdTtr ) );
//         CHECK( params.Tmax == 1200.0 );
//     }

//     SECTION("Testing YAML serialization of ParamsMaierKelleyHKF")
//     {
//         ParamsMaierKelleyHKF params = yaml(params_hkfmk2);
//         CHECK( params.Gf == -39522.064 );
//         CHECK( params.Hf == -31589.2 );
//         CHECK( params.Sr == 143.5112 );
//         CHECK( params.Vr == 3.42e-05 );
//         CHECK( params.a == vec{ 65.39592, 7.610696, 90.3744 } );
//         CHECK( params.b == vec{ 0.0359824, 0.221752, 0.0 } );
//         CHECK( params.c == vec{ 0.0, 0.0, 0.0 } );
//         CHECK( params.Ttr == vec{ 450.0, 620.0 } );
//         CHECK( params.Htr == vec{ 3974.8, 2510.4 } );
//         CHECK( allNaN( params.Vtr ) );
//         CHECK( allNaN( params.dPdTtr ) );
//         CHECK( params.Tmax == 1000.0 );
//     }

//     SECTION("Testing YAML serialization of ParamsMaierKelleyHKF")
//     {
//         ParamsMaierKelleyHKF params = yaml(params_hkfmk3);
//         CHECK( std::isnan( params.Gf ) );
//         CHECK( std::isnan( params.Hf ) );
//         CHECK( params.Sr == 286.604 );
//         CHECK( params.Vr == 0.0001432 );
//         CHECK( params.a == vec{ 369.61456, 409.65544, 488.18912, 461.24416 } );
//         CHECK( params.b == vec{ 0.22643808, 0.14786256, 0.03112896, 0.04217472 } );
//         CHECK( params.c == vec{ -7556304.0, -4167264.0, -1937192.0, -1937192.0 } );
//         CHECK( params.Ttr == vec{ 848.0, 950.0, 1050.0 } );
//         CHECK( allNaN( params.Htr ) );
//         CHECK( allNaN( params.Vtr ) );
//         CHECK( allNaN( params.dPdTtr ) );
//         CHECK( params.Tmax == 1100.0 );
//     }

//     SECTION("Testing YAML serialization of ParamsWaterHKF")
//     {
//         ParamsWaterHKF params = yaml(params_whkf);
//         CHECK( params.Ttr == 273.16 );
//         CHECK( params.Str == 63.312288 );
//         CHECK( params.Gtr == -235517.36 );
//         CHECK( params.Htr == -287721.13 );
//         CHECK( params.Utr == -284039.21 );
//         CHECK( params.Atr == -231856.36 );
//     }

//     SECTION("Testing YAML serialization of ParamsHKF")
//     {
//         ParamsHKF params = yaml(params_hkf);
//         CHECK( params.Gf == 39371.44 );
//         CHECK( params.Hf == -151084.24 );
//         CHECK( params.Sr == 197.4848 );
//         CHECK( params.a1 == 5.8268894e-05 );
//         CHECK( params.a2 == 8251.0572 );
//         CHECK( params.a3 == 0.00049988758 );
//         CHECK( params.a4 == -150377.14 );
//         CHECK( params.c1 == 384.55521 );
//         CHECK( params.c2 == 116047.42 );
//         CHECK( params.wref == -156816.32 );
//         CHECK( params.charge == 0.0 );
//     }

//     SECTION("Testing YAML serialization of ParamsHollandPowell")
//     {
//         ParamsHollandPowell params = yaml(params_hp);
//         CHECK( params.Gf == -4937500.0 );
//         CHECK( params.Hf == -5260650.0 );
//         CHECK( params.Sr == 342.0 );
//         CHECK( params.Vr == 0.00011525 );
//         CHECK( params.a == 677.3 );
//         CHECK( params.b == 0.0 );
//         CHECK( params.c == -3772700.0 );
//         CHECK( params.d == -5044.0 );
//         CHECK( params.alpha == 2.12e-05 );
//         CHECK( params.kappa == 190000000000.0 );
//         CHECK( params.kappap == 2.98 );
//         CHECK( params.kappapp == -1.6e-11 );
//         CHECK( params.numatoms == 20.0 );
//     }

//     SECTION("Testing YAML serialization of ParamsHollandPowellGas")
//     {
//         ParamsHollandPowellGas params = yaml(params_hpg);
//         CHECK( params.Gf == -50710.0 );
//         CHECK( params.Hf == -74810.0 );
//         CHECK( params.Sr == 186.26 );
//         CHECK( params.Vr == 0.0 );
//         CHECK( params.a == 150.1 );
//         CHECK( params.b == 0.002063 );
//         CHECK( params.c == 3427700.0 );
//         CHECK( params.d == -2650.4 );
//     }

//     SECTION("Testing YAML serialization of ParamsHollandPowellLandau")
//     {
//         ParamsHollandPowellLandau params = yaml(params_hpl);
//         CHECK( params.Gf == -2192340.0 );
//         CHECK( params.Hf == -2307040.0 );
//         CHECK( params.Sr == 127.6 );
//         CHECK( params.Vr == 5.16e-05 );
//         CHECK( params.a == 247.5 );
//         CHECK( params.b == -0.003206 );
//         CHECK( params.c == 0.0 );
//         CHECK( params.d == -2051.9 );
//         CHECK( params.alpha == 2.9e-05 );
//         CHECK( params.kappa == 98500000000.0 );
//         CHECK( params.kappap == 4.07 );
//         CHECK( params.kappapp == -4.1e-11 );
//         CHECK( params.numatoms == 7.0 );
//         CHECK( params.Tcr == 1710.0 );
//         CHECK( params.Smax == 10.03 );
//         CHECK( params.Vmax == 5.0e-07 );
//     }
// }
