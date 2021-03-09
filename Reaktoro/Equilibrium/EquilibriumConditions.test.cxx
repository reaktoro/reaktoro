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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumConditions", "[EquilibriumConditions]")
{
    ChemicalSystem system = test::createChemicalSystem();

    EquilibriumSpecs specs(system);

    WHEN("temperature and pressure are input variables - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        EquilibriumConditions conditions(specs);

        conditions.temperature(50, "celsius");
        conditions.pressure(100, "bar");

        auto w = conditions.inputValues();

        CHECK( w.size() == 2 );
        CHECK( w[0] ==  50.0 + 273.15 ); // T in K
        CHECK( w[1] == 100.0 * 1.0e+5 ); // P in Pa

        CHECK_THROWS( conditions.volume(1, "m3") );
    }

    WHEN("temperature and volume are input variables - the Helmholtz energy minimization formulation")
    {
        specs.temperature();
        specs.volume();

        EquilibriumConditions conditions(specs);

        conditions.temperature(40, "celsius");
        conditions.volume(2.0, "m3");

        auto w = conditions.inputValues();

        CHECK( w.size() == 2 );
        CHECK( w[0] == 40.0 + 273.15 ); // T in K
        CHECK( w[1] == 2.0 ); // V in m3

        CHECK_THROWS( conditions.entropy(1, "J/K") );
    }

    WHEN("volume and internal energy are input variables - the entropy maximization formulation")
    {
        specs.volume();
        specs.internalEnergy();

        EquilibriumConditions conditions(specs);

        conditions.volume(1, "cm3");
        conditions.internalEnergy(1, "kJ");

        auto w = conditions.inputValues();

        CHECK( w.size() == 2 );
        CHECK( w[0] == Approx(1.0e-6) ); // V in m3
        CHECK( w[1] == Approx(1.0e+3) ); // U in J

        CHECK_THROWS( conditions.enthalpy(1, "J") );
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

        CHECK( w.size() == 3 );
        CHECK( w[0]  == Approx(35.0 + 273.15) ); // T in K
        CHECK( w[1]  == Approx(23.0 * 1.0e+5) ); // P in Pa
        CHECK( w[2] == 3.5 );                    // pH

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

        CHECK( w.size() == 3 );
        CHECK( w[0] == Approx(2.3 * 1.0e-3) ); // V in m3
        CHECK( w[1] == Approx(1.0e+3) );       // S in J/K
        CHECK( w[2] == Approx(log(40.0)) );    // lnActivity[CO2(g)]

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

        CHECK( w.size() == 6 );
        CHECK( w[0] == Approx(60.0 + 273.15) );  // T in K
        CHECK( w[1] == Approx(200.0 * 1.0e+6) ); // P in Pa
        CHECK( w[2] == Approx(5.0 * 1.0e-9) );   // V in m3
        CHECK( w[3] == Approx(2.0e+6) );         // U in J
        CHECK( w[4] == 2.7 );                    // pH
        CHECK( w[5] == 4.0 );                    // pE

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

        CHECK( w.size() == 4 );
        CHECK( w[0] == 300.0 );      // T in K
        CHECK( w[1] == 1e6 );        // P in Pa
        CHECK( w[2] == K0.value() ); // K0
        CHECK( w[3] == K1.value() ); // K1

        // Check that changing K0 and K1 does not change the input values in object conditions!
        K0 = 7.0;
        K1 = 8.0;

        w = conditions.inputValues();

        CHECK( w[2] == 1.0 ); // K0 at previous value
        CHECK( w[3] == 2.0 ); // K1 at previous value

        // Check that changing input variables "K0" and "K1" in object conditions
        // does not change Param objects K0 and K1!
        conditions.set("K0", 11.0);
        conditions.set("K1", 12.0);

        CHECK( K0.value() == 7.0 ); // recent values of K0 set previously
        CHECK( K1.value() == 8.0 ); // recent values of K1 set previously

        w = conditions.inputValues();

        CHECK( w.size() == 4 );
        CHECK( w[0] == 300.0 ); // corresponding to T
        CHECK( w[1] == 1e6 );   // corresponding to P
        CHECK( w[2] == 11.0 );  // corresponding to K0
        CHECK( w[3] == 12.0 );  // corresponding to K1
    }
}
