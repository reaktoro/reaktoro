// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
using namespace Reaktoro;

TEST_CASE("Testing EquilibriumSolver", "[EquilibriumSolver]")
{
    const auto db = Database({
        Species("H2O"           ).withStandardGibbsEnergy( -237181.72),
        Species("H+"            ).withStandardGibbsEnergy(       0.00),
        Species("OH-"           ).withStandardGibbsEnergy( -157297.48),
        Species("H2"            ).withStandardGibbsEnergy(   17723.42),
        Species("O2"            ).withStandardGibbsEnergy(   16543.54),
        Species("Na+"           ).withStandardGibbsEnergy( -261880.74),
        Species("Cl-"           ).withStandardGibbsEnergy( -131289.74),
        Species("NaCl"          ).withStandardGibbsEnergy( -388735.44),
        Species("HCl"           ).withStandardGibbsEnergy( -127235.44),
        Species("NaOH"          ).withStandardGibbsEnergy( -417981.60),
        Species("Ca++"          ).withStandardGibbsEnergy( -552790.08),
        Species("Mg++"          ).withStandardGibbsEnergy( -453984.92),
        Species("CH4"           ).withStandardGibbsEnergy(  -34451.06),
        Species("CO2"           ).withStandardGibbsEnergy( -385974.00),
        Species("HCO3-"         ).withStandardGibbsEnergy( -586939.89),
        Species("CO3--"         ).withStandardGibbsEnergy( -527983.14),
        Species("CaCl2"         ).withStandardGibbsEnergy( -811696.00),
        Species("CaCO3"         ).withStandardGibbsEnergy(-1099764.40),
        Species("MgCO3"         ).withStandardGibbsEnergy( -998971.84),
        Species("SiO2"          ).withStandardGibbsEnergy( -833410.96),
        Species("CO2(g)"        ).withStandardGibbsEnergy( -394358.74),
        Species("O2(g)"         ).withStandardGibbsEnergy(       0.00),
        Species("H2(g)"         ).withStandardGibbsEnergy(       0.00),
        Species("H2O(g)"        ).withStandardGibbsEnergy( -228131.76),
        Species("CH4(g)"        ).withStandardGibbsEnergy(  -50720.12),
        Species("CO(g)"         ).withStandardGibbsEnergy( -137168.26),
        Species("NaCl(s)"       ).withStandardGibbsEnergy( -384120.49).withName("Halite"   ),
        Species("CaCO3(s)"      ).withStandardGibbsEnergy(-1129177.92).withName("Calcite"  ),
        Species("MgCO3(s)"      ).withStandardGibbsEnergy(-1027833.07).withName("Magnesite"),
        Species("CaMg(CO3)2(s)" ).withStandardGibbsEnergy(-2166307.84).withName("Dolomite" ),
        Species("SiO2(s)"       ).withStandardGibbsEnergy( -856238.86).withName("Quartz"   ),
    });

    // auto aqueous_solution = GENERATE(
    //     AqueousSolution(speciate("H O")),
    //     AqueousSolution(speciate("H O Na Cl")),
    //     AqueousSolution(speciate("H O Na Cl C")),
    //     AqueousSolution(speciate("H O Na Cl C Ca")),
    //     AqueousSolution(speciate("H O Na Cl C Ca Mg")),
    //     AqueousSolution(speciate("H O Na Cl C Ca Mg Si"))
    // );

    // auto gaseous_solution = GENERATE(
    //     GaseousSolution("CO2(g)"),
    //     GaseousSolution("CO2(g) H2O(g)"),
    //     GaseousSolution("CO2(g) H2O(g) O2(g)"),
    //     GaseousSolution("CO2(g) H2O(g) O2(g) H2(g)")
    // );

    // auto minerals = GENERATE(
    //     Minerals("Halite"),
    //     Minerals("Halite Calcite"),
    //     Minerals("Halite Calcite Magnesite"),
    //     Minerals("Halite Calcite Magnesite Dolomite"),
    //     Minerals("Halite Calcite Magnesite Dolomite Quartz")
    // );

    // auto phases = GENERATE_COPY(
    //     Phases(db, aqueous_solution),
    //     Phases(db, aqueous_solution, gaseous_solution),
    //     Phases(db, aqueous_solution, gaseous_solution, minerals),
    //     Phases(db, aqueous_solution, minerals)
    // );

    const auto phases = Phases(db,
        // AqueousSolution(speciate("H O"))
        // AqueousSolution(speciate("H O C"))
        AqueousSolution(speciate("H O Na Cl C Ca Mg Si")),
        GaseousSolution(speciate("H O C"))
        // , Minerals("Halite") // OK
        // , Minerals("Calcite") // OK
        // , Minerals("Magnesite") // OK
        // , Minerals("Dolomite") // OK
        // , Minerals("Quartz") // OK
        // , Minerals("Quartz Calcite") // OK
        // , Minerals("Halite Calcite") // OK
        // , Minerals("Calcite Magnesite") // OK
        // , Minerals("Calcite Magnesite Halite") // ~~FAIL!~~NOW IT WORKS USING zeps to determine stability!
        // , Minerals("Halite Calcite Magnesite Quartz") // WORKS with z > -eps
        // , Minerals("Calcite Dolomite") // OK
        // , Minerals("Magnesite Dolomite") // OK
        // , Minerals("Calcite Magnesite Dolomite") // FAIL!
        // , Minerals("Halite Calcite Magnesite Dolomite")
        , Minerals("Halite Calcite Magnesite Dolomite Quartz")
    );

    const auto system = ChemicalSystem(db, phases);

    // const auto T = GENERATE(25.0, 60.0);
    // const auto P = GENERATE(1.0, 100.0);

    // const auto nH2O    = GENERATE(55.0);
    // const auto nNaCl   = GENERATE(0.0, 1.0, 5.0, 10.0);
    // const auto nCO2    = GENERATE(0.0, 0.1, 1.0, 5.0);
    // const auto nO2     = GENERATE(0.0, 0.1);
    // const auto nH2     = GENERATE(0.0, 0.1);
    // const auto nCH4    = GENERATE(0.0, 0.1);
    // const auto nCaCO3  = GENERATE(0.0, 1e-3, 5.0);
    // const auto nMgCO3  = GENERATE(0.0, 1e-3, 5.0);
    // const auto nSiO2   = GENERATE(0.0, 1e-6, 1.0);

    // const auto T = 60.0;
    // const auto P = 100.0;

    const auto T = 25.0;
    const auto P = 1.0;

    const auto nH2O    = 55.0;
    const auto nNaCl   = 1.0e-2;
    const auto nCO2    = 100.0e-2;
    const auto nO2     = 0.0e-2;
    const auto nH2     = 0.0e-2;
    const auto nCH4    = 0.0e-2;
    const auto nCaCO3  = 1.0e-2;
    const auto nMgCO3  = 2.0e-2;
    const auto nSiO2   = 1.0e-2;

    ChemicalState state(system);
    state.setTemperature(T, "celsius");
    state.setPressure(P, "bar");
    state.setSpeciesAmount("H2O"   , nH2O   , "mol");
    state.setSpeciesAmount("NaCl"  , nNaCl  , "mol");
    state.setSpeciesAmount("CO2"   , nCO2   , "mol");
    // state.setSpeciesAmount("CO2(g)"   , nCO2   , "mol");
    // state.setSpeciesAmount("O2"    , nO2    , "mol");
    // state.setSpeciesAmount("H2"    , nH2    , "mol");
    // state.setSpeciesAmount("CH4"   , nCH4   , "mol");
    state.setSpeciesAmount("CaCO3" , nCaCO3 , "mol");
    // state.setSpeciesAmount("Calcite" , nCaCO3 , "mol");
    state.setSpeciesAmount("MgCO3" , nMgCO3 , "mol");
    state.setSpeciesAmount("SiO2"  , nSiO2  , "mol");

    // VectorXd n(5); n << 55, 2.51011e-13, 2.51011e-13, 2.28356e-31 ,1.14178e-31;
    // state.setSpeciesAmounts(n);


    EquilibriumOptions options;
    // options.optima.output = true;
    options.optima.output = false;
    options.epsilon = 1e-40;
    options.optima.max_iterations = 100;
    options.optima.tolerance = 1e-10;
    // options.optima.tolerance_linear_equality_constraints = options.optima.tolerance;
    // options.optima.tolerance = 1e-20;
    options.optima.kkt.method = Optima::SaddlePointMethod::Rangespace;
    // options.optima.kkt.method = Optima::SaddlePointMethod::Fullspace;
    // options.optima.kkt.method = Optima::SaddlePointMethod::Nullspace;


    options.optima.linesearch.trigger_when_current_error_is_greater_than_initial_error_by_factor = 1;
    options.optima.linesearch.trigger_when_current_error_is_greater_than_previous_error_by_factor = 2;


    EquilibriumSolver solver(system);
    solver.setOptions(options);

    auto start = time();

    ChemicalState copy(state);

    auto result = solver.solve(copy);

    // copy.setSpeciesAmount("CO2", copy.speciesAmount("CO2") + 0.1);

    // solver.solve(copy);

    // copy.setSpeciesAmount("CO2(g)", copy.speciesAmount("CO2(g)") + 0.1);

    // solver.solve(copy);

    // std::cout << "Time (s) = " << elapsed(start) << std::endl;

    CHECK( result.optima.succeeded );
    CHECK( result.optima.iterations <= 48 );
}
