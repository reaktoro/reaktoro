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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPredictor.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
using namespace Reaktoro;

auto createStandardThermoModelH2O() -> StandardThermoModel
{
    return [](real T, real P)
    {
        StandardThermoProps props;
        props.G0 = -237181.72; // in J/mol
        props.V0 = 1.80694885e-5; // in m3/mol
        return props;
    };
}

TEST_CASE("Testing EquilibriumPredictor", "[EquilibriumPredictor]")
{
    const auto db = Database({
        Species("H2O"           ).withStandardThermoModel(createStandardThermoModelH2O()),
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

    Phases phases(db);
    phases.add( AqueousPhase("H2O H+ OH- H2 O2 Na+ Cl- NaCl HCl NaOH") );

    ChemicalSystem system(phases);

    SECTION("when the system is closed, temperature and pressure given")
    {
        EquilibriumSpecs specs(system);
        specs.temperature(); // specify temperature is constrained
        specs.pressure();    // specify pressure is constrained

        EquilibriumConditions conditions0(specs);
        conditions0.temperature(300.0);
        conditions0.pressure(1.0e5);

        ChemicalState state0(system);
        state0.set("H2O" , 55.00, "mol");
        state0.set("NaCl", 0.100, "mol");
        state0.set("O2"  , 0.001, "mol");

        EquilibriumSensitivity sensitivity0(specs);

        EquilibriumSolver solver(specs);
        solver.solve(state0, sensitivity0, conditions0);

        state0.props().update(state0);
        EquilibriumPredictor predictor(state0, sensitivity0);

        ChemicalState state(system);
        state.set("H2O" , 55.50, "mol");
        state.set("NaCl", 0.150, "mol");
        state.set("O2"  , 0.002, "mol");

        const VectorXd naux = state.speciesAmounts();
        const VectorXd b = system.formulaMatrix() * naux;

        EquilibriumConditions conditions(specs);
        conditions.temperature(330.0);
        conditions.pressure(1.1e5);

        predictor.predict(state, conditions);

        const auto dndb0 = sensitivity0.dndb();
        const auto dndw0 = sensitivity0.dndw();
        const auto dpdw0 = sensitivity0.dpdw();
        const auto dqdw0 = sensitivity0.dqdw();
        const auto dpdb0 = sensitivity0.dpdb();
        const auto dqdb0 = sensitivity0.dqdb();
        const auto dudw0 = sensitivity0.dudw();
        const auto dudb0 = sensitivity0.dudb();

        const VectorXd n0 = state0.speciesAmounts();
        const VectorXd p0 = state0.equilibrium().p();
        const VectorXd q0 = state0.equilibrium().q();
        const VectorXd u0 = state0.props();

        const VectorXd w0 = conditions0.inputValues();
        const VectorXd b0 = state0.equilibrium().b();

        const VectorXd w  = conditions.inputValues();

        const VectorXd n = n0 + dndb0*(b - b0) + dndw0*(w - w0);
        const VectorXd p = p0 + dpdb0*(b - b0) + dpdw0*(w - w0);
        const VectorXd q = q0 + dqdb0*(b - b0) + dqdw0*(w - w0);
        const VectorXd u = u0 + dudb0*(b - b0) + dudw0*(w - w0);

        CHECK( n.isApprox(VectorXd(state.speciesAmounts())) );
        CHECK( p.isApprox(VectorXd(state.equilibrium().p())) );
        CHECK( q.isApprox(VectorXd(state.equilibrium().q())) );
        CHECK( u.isApprox(VectorXd(state.props())) );
    }

    SECTION("when the system is closed, temperature and pressure given, O2 is a meta-stable basic species - sensitivity derivatives should be zero")
    {
        EquilibriumSpecs specs(system);
        specs.temperature(); // specify temperature is constrained
        specs.pressure();    // specify pressure is constrained

        EquilibriumConditions conditions0(specs);
        conditions0.temperature(300.0);
        conditions0.pressure(1.0e5);

        ChemicalState state0(system);
        state0.set("H2O" , 55.00, "mol");
        state0.set("NaCl", 0.100, "mol"); // no O2 given here!

        EquilibriumSensitivity sensitivity0(specs);

        EquilibriumSolver solver(specs);
        solver.solve(state0, sensitivity0, conditions0);

        state0.props().update(state0);
        EquilibriumPredictor predictor(state0, sensitivity0);

        ChemicalState state(system);
        state.set("H2O" , 55.50, "mol");
        state.set("NaCl", 0.150, "mol"); // no O2 given here!

        const VectorXd naux = state.speciesAmounts();
        const VectorXd b = system.formulaMatrix() * naux;

        EquilibriumConditions conditions(specs);
        conditions.temperature(330.0);
        conditions.pressure(1.1e5);

        predictor.predict(state, conditions);

        // TODO: Organize tests for EquilibriumPredictor in a way to avoid repeated codes.

        const auto dndb0 = sensitivity0.dndb();
        const auto dndw0 = sensitivity0.dndw();
        const auto dpdw0 = sensitivity0.dpdw();
        const auto dqdw0 = sensitivity0.dqdw();
        const auto dpdb0 = sensitivity0.dpdb();
        const auto dqdb0 = sensitivity0.dqdb();
        const auto dudw0 = sensitivity0.dudw();
        const auto dudb0 = sensitivity0.dudb();

        const VectorXd n0 = state0.speciesAmounts();
        const VectorXd p0 = state0.equilibrium().p();
        const VectorXd q0 = state0.equilibrium().q();
        const VectorXd u0 = state0.props();

        const VectorXd w0 = conditions0.inputValues();
        const VectorXd b0 = state0.equilibrium().b();

        const VectorXd w  = conditions.inputValues();

        const VectorXd n = n0 + dndb0*(b - b0) + dndw0*(w - w0);
        const VectorXd p = p0 + dpdb0*(b - b0) + dpdw0*(w - w0);
        const VectorXd q = q0 + dqdb0*(b - b0) + dqdw0*(w - w0);
        const VectorXd u = u0 + dudb0*(b - b0) + dudw0*(w - w0);

        CHECK( n.isApprox(VectorXd(state.speciesAmounts())) );
        CHECK( p.isApprox(VectorXd(state.equilibrium().p())) );
        CHECK( q.isApprox(VectorXd(state.equilibrium().q())) );
        CHECK( u.isApprox(VectorXd(state.props())) );
    }

    SECTION("when the system is open to H2O and H+, temperature, pressure, volume and pH are given, and basic species O2 is meta-stable")
    {
        EquilibriumSpecs specs(system);
        specs.temperature(); // specify temperature is constrained
        specs.pressure();    // specify pressure is constrained
        specs.volume();      // specify volume is constrained
        specs.pH();          // specify pH is constrained (this creates a q control variable n[H+])
        specs.openTo("H2O"); // specify the system is open to H2O (this creates a p control variable n[H2O])

        EquilibriumConditions conditions0(specs);
        conditions0.temperature(300.0);
        conditions0.pressure(1.0e5);
        conditions0.volume(2.0, "liter"); // H2O will be added as much as needed to fulfil this volume
        conditions0.pH(4.0); // H+ will be added as much as needed to achieve this pH

        ChemicalState state0(system);
        state0.set("H2O" , 55.0, "mol");
        state0.set("NaCl", 0.10, "mol");

        EquilibriumSensitivity sensitivity0(specs);

        EquilibriumSolver solver(specs);
        solver.solve(state0, sensitivity0, conditions0);

        state0.props().update(state0);
        EquilibriumPredictor predictor(state0, sensitivity0);

        ChemicalState state(system);
        state.set("H2O" , 55.5, "mol");
        state.set("NaCl", 0.15, "mol");

        const VectorXd naux = state.speciesAmounts();
        const VectorXd b = system.formulaMatrix() * naux;

        EquilibriumConditions conditions(specs);
        conditions.temperature(330.0);
        conditions.pressure(1.1e5);
        conditions.volume(2.1, "liter"); // H2O will be added as much as needed to fulfil this volume
        conditions.pH(4.2); // H+ will be added as much as needed to achieve this pH

        predictor.predict(state, conditions);

        const auto dndb0 = sensitivity0.dndb();
        const auto dndw0 = sensitivity0.dndw();
        const auto dpdw0 = sensitivity0.dpdw();
        const auto dqdw0 = sensitivity0.dqdw();
        const auto dpdb0 = sensitivity0.dpdb();
        const auto dqdb0 = sensitivity0.dqdb();
        const auto dudw0 = sensitivity0.dudw();
        const auto dudb0 = sensitivity0.dudb();

        const VectorXd n0 = state0.speciesAmounts();
        const VectorXd p0 = state0.equilibrium().p();
        const VectorXd q0 = state0.equilibrium().q();
        const VectorXd u0 = state0.props();

        const VectorXd w0 = conditions0.inputValues();
        const VectorXd b0 = state0.equilibrium().b();

        const VectorXd w  = conditions.inputValues();

        const VectorXd n = n0 + dndb0*(b - b0) + dndw0*(w - w0);
        const VectorXd p = p0 + dpdb0*(b - b0) + dpdw0*(w - w0);
        const VectorXd q = q0 + dqdb0*(b - b0) + dqdw0*(w - w0);
        const VectorXd u = u0 + dudb0*(b - b0) + dudw0*(w - w0);

        CHECK( n.isApprox(VectorXd(state.speciesAmounts())) );
        CHECK( p.isApprox(VectorXd(state.equilibrium().p())) );
        CHECK( q.isApprox(VectorXd(state.equilibrium().q())) );
        CHECK( u.isApprox(VectorXd(state.props())) );
    }
}
