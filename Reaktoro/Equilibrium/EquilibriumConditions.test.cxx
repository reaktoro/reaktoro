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
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumConditions", "[EquilibriumConditions]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto Ns = system.surfaces().size(); // the number of surfaces

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
        conditions.surfaceArea("AqueousPhase:GaseousPhase", 1.0, "m2");
        conditions.surfaceArea("AqueousPhase:Halite", 2.0e6, "mm2");
        conditions.surfaceArea("Calcite", 3.0e4, "cm2");

        CHECK( conditions.inputValue("T")                                      == Approx(298.15) );
        CHECK( conditions.inputValue("P")                                      == Approx(1.0e5) );
        CHECK( conditions.inputValue("pH")                                     == Approx(5.0) );
        CHECK( conditions.inputValue("surfaceArea[AqueousPhase:GaseousPhase]") == Approx(1.0) );
        CHECK( conditions.inputValue("surfaceArea[AqueousPhase:Halite]")       == Approx(2.0) );
        CHECK( conditions.inputValue("surfaceArea[Calcite]")                   == Approx(3.0) );

        conditions.surfaceAreas(ArrayXr{{4.0, 5.0, 6.0}});

        CHECK( conditions.inputValue("surfaceArea[AqueousPhase:GaseousPhase]") == Approx(4.0) );
        CHECK( conditions.inputValue("surfaceArea[AqueousPhase:Halite]")       == Approx(5.0) );
        CHECK( conditions.inputValue("surfaceArea[Calcite]")                   == Approx(6.0) );
    }

    WHEN("temperature and pressure are input variables - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        EquilibriumConditions conditions(specs);

        conditions.temperature(50, "celsius");
        conditions.pressure(100, "bar");

        auto w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == Ns + 2 );

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

        CHECK( conditions.inputValues().size() == Ns + 2 );

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

        CHECK( conditions.inputValues().size() == Ns + 2 );

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

        CHECK( conditions.inputValues().size() == Ns + 3 );

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

        CHECK( conditions.inputValues().size() == Ns + 3 );

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

        CHECK( conditions.inputValues().size() == Ns + 6 );

        CHECK( conditions.inputValue("T")  == Approx(60.0 + 273.15) );  // T in K
        CHECK( conditions.inputValue("P")  == Approx(200.0 * 1.0e+6) ); // P in Pa
        CHECK( conditions.inputValue("V")  == Approx(5.0 * 1.0e-9) );   // V in m3
        CHECK( conditions.inputValue("U")  == Approx(2.0e+6) );         // U in J
        CHECK( conditions.inputValue("pH") == 2.7 );                    // pH
        CHECK( conditions.inputValue("pE") == 4.0 );                    // pE

        CHECK_THROWS( conditions.Eh(14.0, "mV") );
    }

    WHEN("there are Param objects among the input variables")
    {
        // This test exists to ensure that an EquilibriumConditions object does
        // not change the values of the model parameters. These Param objects
        // are wrapper to shared pointers. Changing these objects should be
        // restricted to just a few classes to avoid unexpected side effects.

        Param K0("K0", 1.0);
        Param K1("K1", 2.0);

        specs.temperature();
        specs.pressure();
        specs.addInput(K0);
        specs.addInput(K1);

        EquilibriumConditions conditions(specs);
        conditions.temperature(300.0);
        conditions.pressure(1e6);

        VectorXr w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == Ns + 4 );

        CHECK( conditions.inputValue("T")  == 300.0 );      // T in K
        CHECK( conditions.inputValue("P")  == 1e6 );        // P in Pa
        CHECK( conditions.inputValue("K0") == K0.value() ); // K0
        CHECK( conditions.inputValue("K1") == K1.value() ); // K1

        // Check that changing K0 and K1 does not change the input values in object conditions!
        K0 = 7.0;
        K1 = 8.0;

        w = conditions.inputValues();

        CHECK( conditions.inputValue("K0") == 1.0 ); // K0 at previous value
        CHECK( conditions.inputValue("K1") == 2.0 ); // K1 at previous value

        // Check that changing input variables "K0" and "K1" in object conditions
        // does not change Param objects K0 and K1!
        conditions.set("K0", 11.0);
        conditions.set("K1", 12.0);

        CHECK( K0.value() == 7.0 ); // recent values of K0 set previously
        CHECK( K1.value() == 8.0 ); // recent values of K1 set previously

        w = conditions.inputValues();

        CHECK( conditions.inputValues().size() == Ns + 4 );

        CHECK( conditions.inputValue("T")  == 300.0 ); // corresponding to T
        CHECK( conditions.inputValue("P")  == 1e6 );   // corresponding to P
        CHECK( conditions.inputValue("K0") == 11.0 );  // corresponding to K0
        CHECK( conditions.inputValue("K1") == 12.0 );  // corresponding to K1
    }
}
