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

// C++ includes
#include <iomanip>

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

#define PRINT_INFO_IF_FAILS(x) INFO(#x " = \n" << std::scientific << std::setprecision(16) << x)

// Check if the resulting state of a ChemicalState object in EquilibriumSolver::solve has zero derivative values.
auto checkChemicalEquilibriumStateHasZeroDerivativeValues(ChemicalState const& state)
{
    // Check if temperature, pressure, species amounts have zero derivative seed values
    CHECK(autodiff::grad(state.temperature()) == 0.0);
    CHECK(autodiff::grad(state.pressure()) == 0.0);
    CHECK(autodiff::grad(state.speciesAmounts()).isZero());

    // Check if all chemical and thermodynamic properties have zero derivative seed values
    ArrayStream<real> stream;
    state.props().serialize(stream);

    CHECK(autodiff::grad(ArrayXr(stream)).isZero());
}

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
    options.hessian = GibbsHessian::Exact;
    options.optima.maxiters = 100;
    options.optima.convergence.tolerance = 1e-10;

    EquilibriumResult result;

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

        WHEN("using epsilon 1e-40")
        {
            options.epsilon = 1e-40;
            solver.setOptions(options);

            result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 17 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
        }

        WHEN("using epsilon 1e-16")
        {
            options.epsilon = 1e-16;
            solver.setOptions(options);

            result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 14 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            WHEN("Sensitivity derivatives are considered")
            {
                EquilibriumSensitivity sensitivity;

                result = solver.solve(state, sensitivity); // check a recalculation converges in 0 iterations (even if sensitivity derivatives need to be computed!)

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 0 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

                const auto dndT = sensitivity.dndw("T");
                const auto dndP = sensitivity.dndw("P");
                const auto dndc = sensitivity.dndc();

                PRINT_INFO_IF_FAILS(dndT);
                CHECK( dndT.isApprox(VectorXd({{
                   -2.3437508072185746e-08,
                    2.3437508072185746e-08,
                    2.3437508072185759e-08,
                    0.0000000000000000e+00,
                    0.0000000000000000e+00 }})));

                PRINT_INFO_IF_FAILS(dndP);
                CHECK( dndP.isApprox(VectorXd({{
                    0.0000000000000000e+00,
                    0.0000000000000000e+00,
                    0.0000000000000000e+00,
                    0.0000000000000000e+00,
                    0.0000000000000000e+00 }})));

                PRINT_INFO_IF_FAILS(dndc);
                CHECK( dndc.isApprox(MatrixXd({
                    { 4.9999999384665239e-01,  0.0000000000000000e+00,  1.2306694552322028e-09 },
                    { 6.1533476097769288e-09, -0.0000000000000000e+00,  4.9999999876933054e-01 },
                    { 6.1533476097769296e-09, -0.0000000000000000e+00, -5.0000000123066946e-01 },
                    { 0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00 },
                    { 0.0000000000000000e+00, -0.0000000000000000e+00,  0.0000000000000000e+00 }})));
            }
        }
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

        WHEN("using epsilon 1e-40")
        {
            options.epsilon = 1e-40;
            solver.setOptions(options);

            result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 25 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
        }

        WHEN("using epsilon 1e-16")
        {
            options.epsilon = 1e-16;
            solver.setOptions(options);

            result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 14 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
        }
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

        WHEN("using epsilon 1e-40")
        {
            options.epsilon = 1e-40;
            solver.setOptions(options);

            result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 41 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
        }

        WHEN("using epsilon 1e-16")
        {
            options.epsilon = 1e-16;
            solver.setOptions(options);

            result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 27 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
        }
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

        WHEN("using epsilon 1e-40")
        {
            options.epsilon = 1e-40;
            solver.setOptions(options);

            result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 41 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
        }

        WHEN("using epsilon 1e-16")
        {
            options.epsilon = 1e-16;
            solver.setOptions(options);

            result = solver.solve(state);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 27 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
        }
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

        WHEN("reactivity restrictions are not imposed")
        {
            WHEN("using epsilon 1e-40")
            {
                options.epsilon = 1e-40;
                solver.setOptions(options);

                result = solver.solve(state);

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations <= 62 ); // for some reason, 57 iterations are needed in Linux and 62 in Windows and macOS! This started to appear after ChemicalState was changed (initial species amounts became 1e-16 instead of zero)
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

                result = solver.solve(state); // check a recalculation converges in 0 iterations

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 0 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
            }

            WHEN("using epsilon 1e-16")
            {
                options.epsilon = 1e-16;
                solver.setOptions(options);

                result = solver.solve(state);

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 29 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

                result = solver.solve(state); // check a recalculation converges in 0 iterations

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 0 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);
            }
        }

        WHEN("reactivity restrictions are imposed")
        {
            EquilibriumRestrictions restrictions(system);
            restrictions.cannotIncreaseAbove("Quartz", 0.007, "mol"); // Quartz will precipitate out of 0.01 mol of SiO2(aq) but this will limit to 0.007 mol instead of 0.00973917 mol
            restrictions.cannotDecreaseBelow("MgCO3", 0.10, "mol"); // MgCO3 will be consumed to precipitate Magnesite and Dolomite, but this restriction will prevent it from going below 0.10 moles (without this restriction, it would go to 0.0380553 moles)
            restrictions.cannotReact("Halite"); // the initial amount of Halite, 0.03 mol, would be completely dissolved if this restriction was not imposed

            WHEN("using epsilon 1e-40")
            {
                options.epsilon = 1e-40;
                solver.setOptions(options);

                result = solver.solve(state, restrictions);

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 41 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

                result = solver.solve(state, restrictions); // check a recalculation converges in 0 iterations

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 0 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

                CHECK( state.speciesAmount("Quartz") == Approx(0.007) );
                CHECK( state.speciesAmount("MgCO3")  == Approx(0.100) );
                CHECK( state.speciesAmount("Halite") == Approx(0.030) );
            }

            WHEN("using epsilon 1e-16")
            {
                options.epsilon = 1e-16;
                solver.setOptions(options);

                result = solver.solve(state, restrictions);

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 27 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

                result = solver.solve(state, restrictions); // check a recalculation converges in 0 iterations

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 0 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

                CHECK( state.speciesAmount("Quartz") == Approx(0.007) );
                CHECK( state.speciesAmount("MgCO3")  == Approx(0.100) );
                CHECK( state.speciesAmount("Halite") == Approx(0.030) );
            }
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
        conditions.pH(3.0);

        ChemicalState state(system);
        state.set("H2O", 55, "mol");

        WHEN("using epsilon 1e-40")
        {
            options.epsilon = 1e-40;
            solver.setOptions(options);

            result = solver.solve(state, conditions);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 17 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state, conditions); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            CHECK( state.temperature() == Approx(50.0 + 273.15) );
            CHECK( state.pressure() == Approx(80.0 * 1.0e+5) );
            CHECK( state.speciesAmount("H+") == Approx(0.00099084) );
        }

        WHEN("using epsilon 1e-16")
        {
            options.epsilon = 1e-16;
            solver.setOptions(options);

            result = solver.solve(state, conditions);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 17 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state, conditions); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            CHECK( state.temperature() == Approx(50.0 + 273.15) );
            CHECK( state.pressure() == Approx(80.0 * 1.0e+5) );
            CHECK( state.speciesAmount("H+") == Approx(0.00099084) );

            WHEN("sensitivity derivatives are considered")
            {
                EquilibriumSensitivity sensitivity;

                result = solver.solve(state, sensitivity, conditions); // check a recalculation converges in 0 iterations (even if sensitivity derivatives need to be computed!)

                CHECK( result.optima.succeeded );
                CHECK( result.optima.iterations == 0 );
                checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

                const auto dndT  = sensitivity.dndw("T");
                const auto dndP  = sensitivity.dndw("P");
                const auto dndpH = sensitivity.dndw("pH");
                const auto dndc  = sensitivity.dndc();

                PRINT_INFO_IF_FAILS(dndT);
                CHECK( dndT.isApprox(VectorXd({{
                   -1.1153479330259708e-11,
                   -2.0093291926708909e-16,
                    1.1153479330259708e-11,
                    0.0000000000000000e+00,
                    0.0000000000000000e+00 }})));

                PRINT_INFO_IF_FAILS(dndP);
                CHECK( dndP.isApprox(VectorXd({{
                    0.0000000000000000e+00,
                   -0.0000000000000000e+00,
                   -0.0000000000000000e+00,
                    0.0000000000000000e+00,
                   -0.0000000000000000e+00 }})));

                PRINT_INFO_IF_FAILS(dndpH);
                CHECK( dndpH.isApprox(VectorXd({{
                   -2.7913561385411666e-10,
                   -2.2814928148700963e-03,
                    2.7913561385411666e-10,
                    0.0000000000000000e+00,
                    0.0000000000000000e+00 }})));

                PRINT_INFO_IF_FAILS(dndc);
                CHECK( dndc.isApprox(MatrixXd({
                    { 4.9999999999834693e-01, 0.0000000000000000e+00, -4.9999999999834693e-01 },
                    { 9.0076339999702179e-06, 0.0000000000000000e+00, -9.0076339999702179e-06 },
                    { 1.6530565225577372e-12, 0.0000000000000000e+00, -1.6530565225577372e-12 },
                    { 0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00 },
                    { 0.0000000000000000e+00, 0.0000000000000000e+00,  0.0000000000000000e+00 }})));
            }
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
        conditions.pH(3.0);

        ChemicalState state(system);
        state.set("H2O"  , 55.0, "mol");
        state.set("NaCl" , 0.01, "mol");
        state.set("CO2"  , 10.0, "mol");
        state.set("CaCO3", 0.01, "mol");
        state.set("MgCO3", 0.02, "mol");
        state.set("SiO2" , 0.01, "mol");

        WHEN("using epsilon 1e-40")
        {
            options.epsilon = 1e-40;
            solver.setOptions(options);

            result = solver.solve(state, conditions);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 42 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state, conditions); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            CHECK( state.temperature() == Approx(50.0 + 273.15) );
            CHECK( state.pressure() == Approx(80.0 * 1.0e+5) );
            CHECK( state.speciesAmount("H+") == Approx(0.00099125) );
        }

        WHEN("using epsilon 1e-16")
        {
            options.epsilon = 1e-16;
            solver.setOptions(options);

            result = solver.solve(state, conditions);

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 28 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            result = solver.solve(state, conditions); // check a recalculation converges in 0 iterations

            CHECK( result.optima.succeeded );
            CHECK( result.optima.iterations == 0 );
            checkChemicalEquilibriumStateHasZeroDerivativeValues(state);

            CHECK( state.temperature() == Approx(50.0 + 273.15) );
            CHECK( state.pressure() == Approx(80.0 * 1.0e+5) );
            CHECK( state.speciesAmount("H+") == Approx(0.00099125) );
        }
    }
}
