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
#include <Reaktoro/Common/TimeUtils.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
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

    const auto T = 60.0;  // in celsius
    const auto P = 100.0; // in bar

    EquilibriumOptions options;
    // options.optima.output.active = true;
    options.optima.maxiterations = 100;
    options.optima.convergence.tolerance = 1e-10;

    SECTION("there is only pure water")
    {
        Phases phases(db);
        phases.add( AqueousPhase(speciate("H O")) );

        ChemicalSystem system(phases);

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesAmount("H2O", 55, "mol");

        EquilibriumSolver solver(system);
        solver.setOptions(options);

        auto result = solver.solve(state);

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 15 );

        result = solver.solve(state); // check a recalculation converges in 0 iterations

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 0 );

        SECTION("and sensitivity derivatives are considered")
        {
            EquilibriumSensitivity sensitivity;

            result = solver.solve(state, sensitivity); // check a recalculation converges in 0 iterations (even if sensitivity derivatives need to be computed!)

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );

            const auto dndT = sensitivity.dndw("T");
            const auto dndP = sensitivity.dndw("P");
            const auto dndb = sensitivity.dndb();

            CHECK( dndT == VectorXd({{
                -2.3437508072278416e-08,
                 2.3437508072278416e-08,
                 2.3437508072278413e-08,
                 0.0000000000000000e+00,
                -1.5530798350377789e-34 }}));

            CHECK( dndP == VectorXd({{
                 0.0000000000000000e+00,
                 0.0000000000000000e+00,
                 0.0000000000000000e+00,
                 0.0000000000000000e+00,
                 0.0000000000000000e+00 }}));

            CHECK( dndb == MatrixXd({
                {4.9999999507732190e-01,  0.0000000000000000e+00, -1.1102230246251565e-16},
                {4.9226781169203182e-09, -0.0000000000000000e+00,  5.0000000000000011e-01},
                {4.9226781169203190e-09, -0.0000000000000000e+00, -5.0000000000000000e-01},
                {0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00},
                {4.5454545454545456e-19, -0.0000000000000000e+00, -8.6736173798840353e-35} }));
        }
    }

    SECTION("there is only pure water with allowed extremely tiny species amounts")
    {
        Phases phases(db);
        phases.add( AqueousPhase(speciate("H O")) );

        ChemicalSystem system(phases);

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesAmount("H2O", 55, "mol");

        EquilibriumOptions opts;
        opts = options;
        opts.epsilon = 1e-40;

        EquilibriumSolver solver(system);
        solver.setOptions(opts);

        auto result = solver.solve(state);

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 29 );

        result = solver.solve(state); // check a recalculation converges in 0 iterations

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 0 );
    }

    SECTION("there is only pure water but there are other elements besides H and O with zero amounts")
    {
        Phases phases(db);
        phases.add( AqueousPhase(speciate("H O C Na Cl Ca")) );

        ChemicalSystem system(phases);

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesAmount("H2O", 55, "mol"); // no amount given for species with elements C, Na, Cl, Ca

        EquilibriumSolver solver(system);
        solver.setOptions(options);

        auto result = solver.solve(state);

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 15 );

        result = solver.solve(state); // check a recalculation converges in 0 iterations

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 0 );
    }

    SECTION("there is a more complicated aqueous solution")
    {
        Phases phases(db);
        phases.add( AqueousPhase(speciate("H O Na Cl C Ca Mg Si")) );

        ChemicalSystem system(phases);

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesAmount("H2O"   , 55.0 , "mol");
        state.setSpeciesAmount("NaCl"  , 0.01 , "mol");
        state.setSpeciesAmount("CO2"   , 10.0 , "mol");
        state.setSpeciesAmount("CaCO3" , 0.01 , "mol");
        state.setSpeciesAmount("MgCO3" , 0.02 , "mol");
        state.setSpeciesAmount("SiO2"  , 0.01 , "mol");

        EquilibriumSolver solver(system);
        solver.setOptions(options);

        auto result = solver.solve(state);

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 28 );

        result = solver.solve(state); // check a recalculation converges in 0 iterations

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 0 );
    }

    SECTION("there is an aqueous solution and a gaseous solution")
    {
        Phases phases(db);
        phases.add( AqueousPhase(speciate("H O Na Cl C Ca Mg Si")) );
        phases.add( GaseousPhase(speciate("H O C")) );

        ChemicalSystem system(phases);

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesAmount("H2O"   , 55.0 , "mol");
        state.setSpeciesAmount("NaCl"  , 0.01 , "mol");
        state.setSpeciesAmount("CO2"   , 10.0 , "mol");
        state.setSpeciesAmount("CaCO3" , 0.01 , "mol");
        state.setSpeciesAmount("MgCO3" , 0.02 , "mol");
        state.setSpeciesAmount("SiO2"  , 0.01 , "mol");

        EquilibriumSolver solver(system);
        solver.setOptions(options);

        auto result = solver.solve(state);

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 28 );

        result = solver.solve(state); // check a recalculation converges in 0 iterations

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 0 );
    }

    SECTION("there is an aqueous solution, gaseous solution, several minerals")
    {
        Phases phases(db);
        phases.add( AqueousPhase(speciate("H O Na Cl C Ca Mg Si")) );
        phases.add( GaseousPhase(speciate("H O C")) );
        phases.add( MineralPhases("Halite Calcite Magnesite Dolomite Quartz") );

        ChemicalSystem system(phases);

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesAmount("H2O"   , 55.0 , "mol");
        state.setSpeciesAmount("NaCl"  , 0.01 , "mol");
        state.setSpeciesAmount("CO2"   , 10.0 , "mol");
        state.setSpeciesAmount("CaCO3" , 0.10 , "mol");
        state.setSpeciesAmount("MgCO3" , 0.20 , "mol");
        state.setSpeciesAmount("SiO2"  , 0.01 , "mol");
        state.setSpeciesAmount("Halite", 0.03 , "mol");

        EquilibriumSolver solver(system);
        solver.setOptions(options);

        WHEN("reactivity restrictions are not imposed")
        {
            auto result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 30 );

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
        }

        WHEN("reactivity restrictions are imposed")
        {
            EquilibriumRestrictions restrictions(system);
            restrictions.cannotIncreaseAbove("Quartz", 0.007, "mol"); // Quartz will precipitate out of 0.01 mol of SiO2(aq) but this will limit to 0.007 mol instead of 0.00973917 mol
            restrictions.cannotDecreaseBelow("MgCO3", 0.10, "mol"); // MgCO3 will be consumed to precipitate Magnesite and Dolomite, but this restriction will prevent it from going below 0.10 moles (without this restriction, it would go to 0.0380553 moles)
            restrictions.cannotReact("Halite"); // the initial amount of Halite, 0.03 mol, would be completely dissolved if this restriction was not imposed

            auto result = solver.solve(state, restrictions);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 28 );

            result = solver.solve(state, restrictions); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );

            CHECK( state.speciesAmount("Quartz") == Approx(0.007) );
            CHECK( state.speciesAmount("MgCO3")  == Approx(0.1) );
            CHECK( state.speciesAmount("Halite") == Approx(0.03) );
        }
    }

    SECTION("there is only pure water with given pH")
    {
        Phases phases(db);
        phases.add( AqueousPhase(speciate("H O")) );

        ChemicalSystem system(phases);

        EquilibriumSpecs specs(system);
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumSolver solver(specs);
        solver.setOptions(options);

        EquilibriumConditions conditions(specs);
        conditions.temperature(50.0, "celsius");
        conditions.pressure(80.0, "bar");
        conditions.startWith("H2O", 55, "mol");
        conditions.pH(3.0);

        ChemicalState state(system);

        auto result = solver.solve(state, conditions);

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 17 );

        result = solver.solve(state, conditions); // check a recalculation converges in 0 iterations

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 0 );

        CHECK( state.temperature() == Approx(50.0 + 273.15) );
        CHECK( state.pressure() == Approx(80.0 * 1.0e+5) );
        CHECK( state.speciesAmount("H+") == Approx(0.00099084) );

        SECTION("and sensitivity derivatives are considered")
        {
            EquilibriumSensitivity sensitivity;

            result = solver.solve(state, sensitivity, conditions); // check a recalculation converges in 0 iterations (even if sensitivity derivatives need to be computed!)

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );

            const auto dndT  = sensitivity.dndw("T");
            const auto dndP  = sensitivity.dndw("P");
            const auto dndpH = sensitivity.dndw("pH");
            const auto dndb  = sensitivity.dndb();

            CHECK( dndT == VectorXd({{
                -1.1153479332073403e-11,
                -2.0093291929976317e-16,
                 1.1153479332073403e-11,
                 0.0000000000000000e+00,
                -6.7948459188748201e-35 }}));

            CHECK( dndP == VectorXd({{
                 0.0000000000000000e+00,
                -0.0000000000000000e+00,
                -0.0000000000000000e+00,
                 0.0000000000000000e+00,
                -0.0000000000000000e+00 }}));

            CHECK( dndpH == VectorXd({{
                -2.7913561385420755e-10,
                -2.2814928148700941e-03,
                 2.7913561385420755e-10,
                 0.0000000000000000e+00,
                 0.0000000000000000e+00 }}));

            CHECK( dndb == MatrixXd({
                { 4.9999999999889794e-01,  0.0000000000000000e+00, -4.9999999999889794e-01 },
                { 9.0076339999801367e-06, -0.0000000000000000e+00, -9.0076339999801367e-06 },
                { 1.1020442994623832e-12, -0.0000000000000000e+00, -1.1020442994623832e-12 },
                { 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00 },
                { 4.5454545454545456e-19, -0.0000000000000000e+00, -4.5454545454545456e-19 }}));
        }
    }

    SECTION("there is an aqueous solution with given pH in equilibrium with a gaseous solution")
    {
        Phases phases(db);
        phases.add( AqueousPhase(speciate("H O Na Cl C Ca Mg Si")) );
        phases.add( GaseousPhase(speciate("H O C")) );

        ChemicalSystem system(phases);

        EquilibriumSpecs specs(system);
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumSolver solver(specs);
        solver.setOptions(options);

        EquilibriumConditions conditions(specs);
        conditions.temperature(50.0, "celsius");
        conditions.pressure(80.0, "bar");
        conditions.startWith("H2O"   , 55.0 , "mol");
        conditions.startWith("NaCl"  , 0.01 , "mol");
        conditions.startWith("CO2"   , 10.0 , "mol");
        conditions.startWith("CaCO3" , 0.01 , "mol");
        conditions.startWith("MgCO3" , 0.02 , "mol");
        conditions.startWith("SiO2"  , 0.01 , "mol");
        conditions.pH(3.0);

        ChemicalState state(system);

        auto result = solver.solve(state, conditions);

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 28 );

        result = solver.solve(state, conditions); // check a recalculation converges in 0 iterations

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 0 );

        CHECK( state.temperature() == Approx(50.0 + 273.15) );
        CHECK( state.pressure() == Approx(80.0 * 1.0e+5) );
        CHECK( state.speciesAmount("H+") == Approx(0.00099125) );
    }
}
