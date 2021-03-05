// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPredictor.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
using namespace Reaktoro;

TEST_CASE("Testing EquilibriumPredictor", "[EquilibriumPredictor]")
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

    Phases phases(db);
    phases.add( AqueousPhase("H2O H+ OH- H2 O2") );

    ChemicalSystem system(phases);

    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();

    EquilibriumConditions conditions0(specs);
    conditions0.temperature(300.0);
    conditions0.pressure(1.0e5);
    conditions0.startWith("H2O", 55.0, "mol");

    ChemicalState state0(system);
    EquilibriumSensitivity sensitivity0(specs);

    EquilibriumSolver solver(specs);
    solver.solve(state0, sensitivity0, conditions0);

    EquilibriumPredictor predictor(state0, sensitivity0);

    ChemicalState state(system);

    EquilibriumConditions conditions(specs);
    conditions.temperature(330.0);
    conditions.pressure(1.1e5);
    conditions.startWith("H2O", 55.1, "mol");

    predictor.predict(state, conditions);

    const auto dndb = sensitivity0.dndb();
    const auto dndw = sensitivity0.dndw();
    const auto dpdw = sensitivity0.dpdw();
    const auto dqdw = sensitivity0.dqdw();
    const auto dpdb = sensitivity0.dpdb();
    const auto dqdb = sensitivity0.dqdb();
    const auto dudw = sensitivity0.dudw();
    const auto dudb = sensitivity0.dudb();

    // conditions.initialComponentAmountsCompute(n)
    // conditions.initialComponentAmountsComputeOrRetrieve(n)

    // EquilibriumDims dims(specs);

    // const MatrixXd dndb = random(dims.Nn, dims.Nb);
    // const MatrixXd dndw = random(dims.Nn, dims.Nw);
    // const MatrixXd dpdw = random(dims.Np, dims.Nw);
    // const MatrixXd dqdw = random(dims.Nq, dims.Nw);
    // const MatrixXd dndb = random(dims.Nn, dims.Nb);
    // const MatrixXd dpdb = random(dims.Np, dims.Nb);
    // const MatrixXd dqdb = random(dims.Nq, dims.Nb);
    // const MatrixXd dudw = random(dims.Nu, dims.Nw);
    // const MatrixXd dudb = random(dims.Nu, dims.Nb);

    // EquilibriumSensitivity sensitivity0(specs);
    // sensitivity0.dndb(dndb);
    // sensitivity0.dndw(dndw);
    // sensitivity0.dpdw(dpdw);
    // sensitivity0.dqdw(dqdw);
    // sensitivity0.dndb(dndb);
    // sensitivity0.dpdb(dpdb);
    // sensitivity0.dqdb(dqdb);
    // sensitivity0.dudw(dudw);
    // sensitivity0.dudb(dudb);

    // EquilibriumPredictor predictor(state0, sensitivity0);


}
