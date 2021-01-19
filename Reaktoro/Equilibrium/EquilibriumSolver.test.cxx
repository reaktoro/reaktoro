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

    const auto T = 60.0;  // in celsius
    const auto P = 100.0; // in bar

    EquilibriumOptions options;
    // options.optima.output.active = true;
    options.optima.maxiterations = 100;
    options.optima.convergence.tolerance = 1e-10;

    WHEN("there is only pure water")
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
        CHECK( result.optima.iterations == 14 );
    }

    WHEN("there is only pure water with allowed extremely tiny species amounts")
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
        CHECK( result.optima.iterations == 28 );
    }

    WHEN("there is only pure water but there are other elements besides H and O with zero amounts")
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
        CHECK( result.optima.iterations == 14 );
    }

    WHEN("there is a more complicated aqueous solution")
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
        CHECK( result.optima.iterations == 27 );
    }

    WHEN("there is an aqueous solution and a gaseous solution")
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
    }

    WHEN("there is an aqueous solution, gaseous solution, several minerals")
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
        state.setSpeciesAmount("CaCO3" , 0.01 , "mol");
        state.setSpeciesAmount("MgCO3" , 0.02 , "mol");
        state.setSpeciesAmount("SiO2"  , 0.01 , "mol");

        EquilibriumSolver solver(system);
        solver.setOptions(options);

        auto result = solver.solve(state);

        CHECK( result.optima.succeeded );
        CHECK( result.optima.iterations == 28 );
    }
}
