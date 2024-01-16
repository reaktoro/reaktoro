// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumConditions", "[EquilibriumConditions]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto Nn = system.species().size();  // the number of species
    const auto Ne = system.elements().size(); // the number of elements
    const auto Nb = Ne + 1;                   // the number of elements and charge

    EquilibriumSpecs specs(system);

    WHEN("the EquilibriumConditions object holds a default state")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumConditions conditions(specs);

        // Check all input values are initially set to NaN
        for(auto value : conditions.inputValues() )
            CHECK( std::isnan(value.val()) );

        conditions.temperature(25.0, "celsius");
        conditions.pressure(1.0, "bar");
        conditions.pH(5.0);

        CHECK( conditions.inputValue("T")  == Approx(298.15) );
        CHECK( conditions.inputValue("P")  == Approx(1.0e5) );
        CHECK( conditions.inputValue("pH") == Approx(5.0) );
    }

    WHEN("the EquilibriumConditions object holds a default state and inputValuesGetOrCompute is used")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumConditions conditions(specs);

        ChemicalState state0(system);
        conditions.temperature(40.0, "celsius");
        conditions.pressure(13.0, "bar");

        CHECK_THROWS( conditions.inputValuesGetOrCompute(state0) ); // pH has not yet been specified, and it cannot be fetched from state0

        conditions.pH(9.0);

        const auto w = conditions.inputValuesGetOrCompute(state0);

        conditions.setInputVariables(w);

        CHECK( conditions.inputValue("T")  == Approx(273.15 + 40.0) );
        CHECK( conditions.inputValue("P")  == Approx(13.0e5) );
        CHECK( conditions.inputValue("pH") == Approx(9.0) );
    }

    WHEN("temperature and pressure are input variables - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        EquilibriumConditions conditions(specs);

        conditions.temperature(50, "celsius");
        conditions.pressure(100, "bar");

        auto w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == 2 );

        CHECK( conditions.inputValue("T") ==  50.0 + 273.15 ); // T in K
        CHECK( conditions.inputValue("P") == 100.0 * 1.0e+5 ); // P in Pa

        CHECK_THROWS( conditions.volume(1, "m3") );

        CHECK( conditions.lowerBoundsControlVariablesP().size() == 0 ); // there are no p control variables
        CHECK( conditions.upperBoundsControlVariablesP().size() == 0 ); // there are no p control variables
    }

    WHEN("temperature and volume are input variables - the Helmholtz energy minimization formulation")
    {
        specs.temperature();
        specs.volume();

        EquilibriumConditions conditions(specs);

        conditions.temperature(40, "celsius");
        conditions.volume(2.0, "m3");

        auto w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == 2 );

        CHECK( conditions.inputValue("T") == 40.0 + 273.15 ); // T in K
        CHECK( conditions.inputValue("V") == 2.0 ); // V in m3

        CHECK_THROWS( conditions.entropy(1, "J/K") );

        const auto plower = conditions.lowerBoundsControlVariablesP();
        const auto pupper = conditions.upperBoundsControlVariablesP();

        REQUIRE( plower.size() == 1 ); // there is only one p control variable: P
        REQUIRE( pupper.size() == 1 ); // there is only one p control variable: P

        CHECK(( std::isinf(plower[0]) && plower[0] < 0.0 )); // default lower bound for P: -inf
        CHECK(( std::isinf(pupper[0]) && pupper[0] > 0.0 )); // default upper bound for P: inf

        conditions.setLowerBoundPressure(0.01, "bar");
        conditions.setUpperBoundPressure(1.00, "GPa");

        CHECK( plower[0] == 0.01e+5 ); // check new lower bound for P
        CHECK( pupper[0] == 1.00e+9 ); // check new upper bound for P
    }

    WHEN("volume and internal energy are input variables - the entropy maximization formulation")
    {
        specs.volume();
        specs.internalEnergy();

        EquilibriumConditions conditions(specs);

        conditions.volume(1, "cm3");
        conditions.internalEnergy(1, "kJ");

        auto w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == 2 );

        CHECK( conditions.inputValue("V") == Approx(1.0e-6) ); // V in m3
        CHECK( conditions.inputValue("U") == Approx(1.0e+3) ); // U in J

        CHECK_THROWS( conditions.enthalpy(1, "J") );

        const auto plower = conditions.lowerBoundsControlVariablesP();
        const auto pupper = conditions.upperBoundsControlVariablesP();

        REQUIRE( plower.size() == 2 ); // there are two p control variables: T and P
        REQUIRE( pupper.size() == 2 ); // there are two p control variables: T and P

        CHECK(( std::isinf(plower[0]) && plower[0] < 0.0 )); // default lower bound for T: -inf
        CHECK(( std::isinf(plower[1]) && plower[1] < 0.0 )); // default lower bound for P: -inf

        CHECK(( std::isinf(pupper[0]) && pupper[0] > 0.0 )); // default upper bound for T: inf
        CHECK(( std::isinf(pupper[1]) && pupper[1] > 0.0 )); // default upper bound for P: inf

        conditions.setLowerBoundTemperature(0.000, "celsius");
        conditions.setUpperBoundTemperature(100.0, "celsius");

        CHECK( plower[0] == Approx(0.000 + 273.15) ); // check new lower bound for T
        CHECK( pupper[0] == Approx(100.0 + 273.15) ); // check new upper bound for T

        conditions.setLowerBoundPressure(0.01, "bar");
        conditions.setUpperBoundPressure(1.00, "GPa");

        CHECK( plower[1] == Approx(0.01e+5) ); // check new lower bound for P
        CHECK( pupper[1] == Approx(1.00e+9) ); // check new upper bound for P
    }

    WHEN("temperature, pressure, and pH are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumConditions conditions(specs);

        conditions.temperature(35, "celsius");
        conditions.pressure(23, "bar");
        conditions.pH(3.5);

        auto w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == 3 );

        CHECK( conditions.inputValue("T")  == Approx(35.0 + 273.15) ); // T in K
        CHECK( conditions.inputValue("P")  == Approx(23.0 * 1.0e+5) ); // P in Pa
        CHECK( conditions.inputValue("pH") == 3.5 );                   // pH

        CHECK_THROWS( conditions.pE(5.0) );
    }

    WHEN("volume, entropy, and activity[CO2(g)] are input variables")
    {
        specs.volume();
        specs.entropy();
        specs.activity("CO2(g)");

        EquilibriumConditions conditions(specs);

        conditions.volume(2.3, "dm3");
        conditions.entropy(1.0, "kJ/K");
        conditions.activity("CO2(g)", 40.0);

        auto w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == 3 );

        CHECK( conditions.inputValue("V")             == Approx(2.3 * 1.0e-3) ); // V in m3
        CHECK( conditions.inputValue("S")             == Approx(1.0e+3) );       // S in J/K
        CHECK( conditions.inputValue("ln(a[CO2(g)])") == Approx(log(40.0)) );    // ln(a[CO2(g)])

        CHECK_THROWS( conditions.chemicalPotential("H2O(aq)", 100.0, "J/mol") );
    }

    WHEN("temperature, pressure, volume, internal energy, pH, and pE are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.volume();
        specs.internalEnergy();
        specs.pH();
        specs.pE();
        specs.openTo("CO2");
        specs.openTo("CH4");

        EquilibriumConditions conditions(specs);

        conditions.temperature(60, "celsius");
        conditions.pressure(200, "MPa");
        conditions.volume(5.0, "mm3");
        conditions.internalEnergy(2.0, "MJ");
        conditions.pH(2.7);
        conditions.pE(4.0);

        auto w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == 6 );

        CHECK( conditions.inputValue("T")  == Approx(60.0 + 273.15) );  // T in K
        CHECK( conditions.inputValue("P")  == Approx(200.0 * 1.0e+6) ); // P in Pa
        CHECK( conditions.inputValue("V")  == Approx(5.0 * 1.0e-9) );   // V in m3
        CHECK( conditions.inputValue("U")  == Approx(2.0e+6) );         // U in J
        CHECK( conditions.inputValue("pH") == 2.7 );                    // pH
        CHECK( conditions.inputValue("pE") == 4.0 );                    // pE

        CHECK_THROWS( conditions.Eh(14.0, "mV") );
    }

    WHEN("there are custom inputs among the input variables")
    {
        real K0 = 1.0;
        real K1 = 2.0;

        specs.temperature();
        specs.pressure();
        specs.addInput("K0");
        specs.addInput("K1");

        EquilibriumConditions conditions(specs);
        conditions.temperature(300.0);
        conditions.pressure(1e6);
        conditions.set("K0", K0);
        conditions.set("K1", K1);

        VectorXr w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == 4 );

        CHECK( conditions.inputValue("T")  == 300.0 ); // T in K
        CHECK( conditions.inputValue("P")  == 1e6 );   // P in Pa
        CHECK( conditions.inputValue("K0") == K0 );    // K0
        CHECK( conditions.inputValue("K1") == K1 );    // K1
    }

    WHEN("the initial amounts of conservative components are set")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();
        specs.addReactivityConstraint({ "Reaction1", VectorXd::Random(Nn).cwiseAbs(), {} });
        specs.addReactivityConstraint({ "Reaction2", VectorXd::Random(Nn).cwiseAbs(), {} });

        EquilibriumConditions conditions(specs);

        CHECK( conditions.initialComponentAmounts().size() == 0 );

        WHEN("using setInitialComponentAmounts")
        {
            const ArrayXd b0  = ArrayXd::Random(Nb).cwiseAbs();
            const ArrayXd xi0 = ArrayXd::Random(2).cwiseAbs(); // two reactivity constraints

            ArrayXd c0(Nb + 2);
            c0 << b0, xi0;

            CHECK_NOTHROW( conditions.setInitialComponentAmounts(c0) );
            CHECK( conditions.initialComponentAmounts().isApprox(c0) );
        }

        WHEN("using setInitialComponentAmountsFromSpeciesAmounts")
        {
            const ArrayXd n0 = ArrayXd::Random(Nn).cwiseAbs();
            const ArrayXd c0 = conditions.initialComponentAmountsGetOrCompute(n0);

            CHECK_NOTHROW( conditions.setInitialComponentAmountsFromSpeciesAmounts(n0) );
            CHECK( conditions.initialComponentAmounts().isApprox(c0) );
        }

        WHEN("using setInitialComponentAmountsFromState")
        {
            const ArrayXd n0 = ArrayXd::Random(Nn).cwiseAbs();
            const ArrayXd c0 = conditions.initialComponentAmountsGetOrCompute(n0);

            ChemicalState state0(system);
            state0.setSpeciesAmounts(n0);

            CHECK_NOTHROW( conditions.setInitialComponentAmountsFromState(state0) );
            CHECK( conditions.initialComponentAmounts().isApprox(c0) );
        }
    }
}
