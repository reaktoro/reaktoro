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

// ThermoFun includes
#include <ThermoFun/Database.h>

// Reaktoro includes
#include <Reaktoro/Extensions/ThermoFun/ThermoFunDatabase.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ThermoFunDatabase", "[ThermoFunDatabase]")
{
    ThermoFunDatabase::disableLogging(); // Disable logs from ThermoFun for the tests below that are expected to overwrite species or raise exceptions

    const auto T = 298.15;
    const auto P = 1.0e5;

    auto const aq17path = REAKTORO_DATABASES_DIR"/thermofun/aq17-thermofun.json";
    auto const cemdata18path = REAKTORO_DATABASES_DIR"/thermofun/cemdata18-thermofun.json";

    Species species;
    StandardThermoProps props;

    //--------------------------------------------------------------------------------------------------
    // Testing constructor ThermoFunDatabase(name)
    //--------------------------------------------------------------------------------------------------
    CHECK_NOTHROW( ThermoFunDatabase("aq17") );
    CHECK_NOTHROW( ThermoFunDatabase("cemdata18") );
    CHECK_NOTHROW( ThermoFunDatabase("heracles") );
    CHECK_NOTHROW( ThermoFunDatabase("mines16") );
    CHECK_NOTHROW( ThermoFunDatabase("psinagra-12-07") );
    CHECK_NOTHROW( ThermoFunDatabase("slop98-organic") );
    CHECK_NOTHROW( ThermoFunDatabase("slop98") );

    //--------------------------------------------------------------------------------------------------
    // Testing construction of ThermoFunDatabase with other methods
    //--------------------------------------------------------------------------------------------------
    WHEN("Constructing ThermoFunDatabase using ThermoFunDatabase(filename)::fromFiles")
    {
        auto const db1 = ThermoFunDatabase("aq17");
        auto const db2 = ThermoFunDatabase(ThermoFun::Database(aq17path));
        auto const db3 = ThermoFunDatabase::fromFile(aq17path);
        auto const db4 = ThermoFunDatabase::fromFiles({ aq17path });

        Vec<ThermoFunDatabase> dbs = {db1, db2, db3, db4};

        for(auto const& db : dbs)
        {
            //--------------------------------------------------------------------------------------------------
            // Testing attributes and thermodynamic properties of H2O@
            //--------------------------------------------------------------------------------------------------
            species = db.species().get("H2O@");
            CHECK( species.formula().equivalent("H2O") );
            CHECK( species.substance() == "Water HGK" );
            CHECK( species.aggregateState() == AggregateState::Aqueous );
            CHECK( species.charge() == 0 );
            CHECK( species.molarMass() == Approx(0.0180153) );

            props = species.standardThermoProps(T, P);
            CHECK( props.G0  == Approx(-2.371817e+05) );
            CHECK( props.H0  == Approx(-2.858310e+05) );
            CHECK( props.V0  == Approx( 1.806862e-05) );
            CHECK( props.VT0 == Approx(0.00000000000) );
            CHECK( props.VP0 == Approx(0.00000000000) );
            CHECK( props.Cp0 == Approx( 7.532758e+01) );

            //--------------------------------------------------------------------------------------------------
            // Testing attributes and thermodynamic properties of CO3-2
            //--------------------------------------------------------------------------------------------------
            species = db.species().get("CO3-2");
            CHECK( species.formula().equivalent("CO3-2") );
            CHECK( species.substance() == "CO3-2 carbonate ion" );
            CHECK( species.aggregateState() == AggregateState::Aqueous );
            CHECK( species.charge() == -2 );
            CHECK( species.molarMass() == Approx(0.0600100979) );

            props = species.standardThermoProps(T, P);
            CHECK( props.G0  == Approx(-5.279830e+05) );
            CHECK( props.H0  == Approx(-6.752359e+05) );
            CHECK( props.V0  == Approx(-6.063738e-06) );
            CHECK( props.VT0 == Approx(0.00000000000) );
            CHECK( props.VP0 == Approx(0.00000000000) );
            CHECK( props.Cp0 == Approx(-3.228612e+02) );

            //--------------------------------------------------------------------------------------------------
            // Testing attributes and thermodynamic properties of Ca+2
            //--------------------------------------------------------------------------------------------------
            species = db.species().get("Ca+2");
            CHECK( species.formula().equivalent("Ca+2") );
            CHECK( species.substance() == "Ca+2 ion" );
            CHECK( species.aggregateState() == AggregateState::Aqueous );
            CHECK( species.charge() == +2 );
            CHECK( species.molarMass() == Approx(0.040076902) );

            props = species.standardThermoProps(T, P);
            CHECK( props.G0  == Approx(-5.528210e+05) );
            CHECK( props.H0  == Approx(-5.431003e+05) );
            CHECK( props.V0  == Approx(-1.844093e-05) );
            CHECK( props.VT0 == Approx(0.00000000000) );
            CHECK( props.VP0 == Approx(0.00000000000) );
            CHECK( props.Cp0 == Approx(-3.099935e+01) );

            //--------------------------------------------------------------------------------------------------
            // Testing attributes and thermodynamic properties of CO2
            //--------------------------------------------------------------------------------------------------
            species = db.species().get("CO2");
            CHECK( species.formula().equivalent("CO2") );
            CHECK( species.substance() == "Carbon dioxide (CO2)" );
            CHECK( species.aggregateState() == AggregateState::Gas );
            CHECK( species.charge() == 0 );
            CHECK( species.molarMass() == Approx(0.0440096006) );

            props = species.standardThermoProps(T, P);
            CHECK( props.G0  == Approx(-3.943510e+05) );
            CHECK( props.H0  == Approx(-3.935472e+05) );
            CHECK( props.V0  == Approx(0.00000000000) );
            CHECK( props.VT0 == Approx(0.00000000000) );
            CHECK( props.VP0 == Approx(0.00000000000) );
            CHECK( props.Cp0 == Approx( 3.710812e+01) );

            //--------------------------------------------------------------------------------------------------
            // Testing attributes and thermodynamic properties of Calcite
            //--------------------------------------------------------------------------------------------------
            species = db.species().get("Calcite");
            CHECK( species.formula().equivalent("CaCO3") );
            CHECK( species.substance() == "Calcite (cc)" );
            CHECK( species.aggregateState() == AggregateState::CrystallineSolid );
            CHECK( species.charge() == 0 );
            CHECK( species.molarMass() == Approx(0.1000869999) );

            props = species.standardThermoProps(T, P);
            CHECK( props.G0  == Approx(-1.129195e+06) );
            CHECK( props.H0  == Approx(-1.207470e+06) );
            CHECK( props.V0  == Approx( 3.689000e-05) );
            CHECK( props.VT0 == Approx(0.00000000000) );
            CHECK( props.VP0 == Approx(0.00000000000) );
            CHECK( props.Cp0 == Approx( 8.337073e+01) );
        }
    }

    WHEN("Constructing ThermoFunDatabase using ThermoFunDatabase::fromFiles")
    {
        auto const db = ThermoFunDatabase::fromFiles({aq17path, cemdata18path});

        //--------------------------------------------------------------------------------------------------
        // Testing attributes and thermodynamic properties of H2O@, which was replaced by
        //--------------------------------------------------------------------------------------------------
        species = db.species().get("H2O@");
        CHECK( species.formula().equivalent("H2O") );
        CHECK( species.substance() == "H2O  l" ); // instead of Water HGK
        CHECK( species.aggregateState() == AggregateState::Aqueous );
        CHECK( species.charge() == 0 );
        CHECK( species.molarMass() == Approx(0.0180153) );

        props = species.standardThermoProps(T, P);
        CHECK( props.G0  == Approx(-2.371817e+05) );
        CHECK( props.H0  == Approx(-2.858310e+05) );
        CHECK( props.V0  == Approx( 1.806862e-05) );
        CHECK( props.VT0 == Approx(0.00000000000) );
        CHECK( props.VP0 == Approx(0.00000000000) );
        CHECK( props.Cp0 == Approx( 7.532758e+01) );

        //--------------------------------------------------------------------------------------------------
        // Testing attributes and thermodynamic properties of Fe-hemicarbonate
        //   Note: Applicable when cemdata18 was attached with ThermoFunDatabase::fromFiles
        //--------------------------------------------------------------------------------------------------

        species = db.species().get("Fe-hemicarbonate");
        CHECK( species.formula().equivalent("Ca3O3Fe2O3(CaCO3)0.5(CaO2H2)0.5(H2O)9.5") );
        CHECK( species.substance() == "Fe-hemicarbonate" );
        CHECK( species.aggregateState() == AggregateState::CrystallineSolid );
        CHECK( species.charge() == 0 );
        CHECK( species.molarMass() == Approx(0.5861556005) );

        props = species.standardThermoProps(T, P);
        CHECK( props.G0  == Approx(-5952868.4) );
        CHECK( props.H0  == Approx(-6581020.3408304) );
        CHECK( props.V0  == Approx(27.33930015564e-05) );
        CHECK( props.VT0 == Approx(0.00000000000) );
        CHECK( props.VP0 == Approx(0.00000000000) );
        CHECK( props.Cp0 == Approx(841.18804931641) );
    }

    CHECK_THROWS( ThermoFunDatabase("not-a-valid-file-name") );
    CHECK_THROWS( ThermoFunDatabase::withName("not-a-valid-file-name") );
    CHECK_THROWS( ThermoFunDatabase::fromFile("not-a-valid-file-name") );
}
