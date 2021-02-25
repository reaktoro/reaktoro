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

    WHEN("temperature and pressure are input parameters - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        EquilibriumConditions conditions(specs);

        conditions.temperature(50, "celsius");
        conditions.pressure(100, "bar");

        auto params = conditions.params();

        CHECK( params.size() == 2 );
        CHECK( params.get("T").value() ==  50.0 + 273.15 ); // in K
        CHECK( params.get("P").value() == 100.0 * 1.0e+5 ); // in Pa

        CHECK_THROWS( conditions.volume(1, "m3") );
    }

    WHEN("temperature and volume are input parameters - the Helmholtz energy minimization formulation")
    {
        specs.temperature();
        specs.volume();

        EquilibriumConditions conditions(specs);

        conditions.temperature(40, "celsius");
        conditions.volume(2.0, "m3");

        auto params = conditions.params();

        CHECK( params.size() == 2 );
        CHECK( params.get("T").value() == 40.0 + 273.15 ); // in K
        CHECK( params.get("V").value() == 2.0 ); // in m3

        CHECK_THROWS( conditions.entropy(1, "J/K") );
    }

    WHEN("volume and internal energy are input parameters - the entropy maximization formulation")
    {
        specs.volume();
        specs.internalEnergy();

        EquilibriumConditions conditions(specs);

        conditions.volume(1, "cm3");
        conditions.internalEnergy(1, "kJ");

        auto params = conditions.params();

        CHECK( params.size() == 2 );
        CHECK( params.get("V").value() == Approx(1.0e-6) ); // in m3
        CHECK( params.get("U").value() == Approx(1.0e+3) ); // in J

        CHECK_THROWS( conditions.enthalpy(1, "J") );
    }

    WHEN("temperature, pressure, and pH are input parameters")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        EquilibriumConditions conditions(specs);

        conditions.temperature(35, "celsius");
        conditions.pressure(23, "bar");
        conditions.pH(3.5);

        auto params = conditions.params();

        CHECK( params.size() == 3 );
        CHECK( params.get("T").value()  == Approx(35.0 + 273.15) ); // in K
        CHECK( params.get("P").value()  == Approx(23.0 * 1.0e+5) ); // in Pa
        CHECK( params.get("pH").value() == 3.5 );

        CHECK_THROWS( conditions.pE(5.0) );
    }

    WHEN("volume, entropy, and activity[CO2(g)] are input parameters")
    {
        specs.volume();
        specs.entropy();
        specs.activity("CO2(g)");

        EquilibriumConditions conditions(specs);

        conditions.volume(2.3, "dm3");
        conditions.entropy(1.0, "kJ/K");
        conditions.activity("CO2(g)", 40.0);

        auto params = conditions.params();

        CHECK( params.size() == 3 );
        CHECK( params.get("V").value() == Approx(2.3 * 1.0e-3) ); // in m3
        CHECK( params.get("S").value() == Approx(1.0e+3) );       // in J/K
        CHECK( params.get("lnActivity[CO2(g)]").value() == Approx(log(40.0)) );

        CHECK_THROWS( conditions.chemicalPotential("H2O(aq)", 100.0, "J/mol") );
    }

    WHEN("temperature, pressure, volume, internal energy, pH, and pE are input parameters")
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

        auto params = conditions.params();

        CHECK( params.size() == 6 );
        CHECK( params.get("T").value()  == Approx(60.0 + 273.15) );  // in K
        CHECK( params.get("P").value()  == Approx(200.0 * 1.0e+6) ); // in Pa
        CHECK( params.get("V").value()  == Approx(5.0 * 1.0e-9) );   // in m3
        CHECK( params.get("U").value()  == Approx(2.0e+6) );         // in J
        CHECK( params.get("pH").value() == 2.7 );
        CHECK( params.get("pE").value() == 4.0 );

        CHECK_THROWS( conditions.Eh(14.0, "mV") );
    }
}
