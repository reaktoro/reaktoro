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
#include <Reaktoro/Core/Database.hpp>
using namespace Reaktoro;

namespace test {

/// Return mock standard thermodynamic model for an aqueous species.
auto standardThermoModelAqueous(real T, real P) -> StandardThermoProps
{
    StandardThermoProps props;
    props.G0  = 2.1 * (T*P)*(T*P);
    props.H0  = 2.2 * (T*P)*(T*P);
    props.V0  = 2.3 * (T*P)*(T*P);
    props.Cp0 = 2.4 * (T*P)*(T*P);
    props.Cv0 = 2.5 * (T*P)*(T*P);
    return props;
};

/// Return mock standard thermodynamic model for a gaseous species.
auto standardThermoModelGaseous(real T, real P) -> StandardThermoProps
{
    StandardThermoProps props;
    props.G0  = 0.1 * (T*P)*(T*P);
    props.H0  = 0.2 * (T*P)*(T*P);
    props.V0  = 0.3 * (T*P)*(T*P);
    props.Cp0 = 0.4 * (T*P)*(T*P);
    props.Cv0 = 0.5 * (T*P)*(T*P);
    return props;
};

/// Return mock standard thermodynamic model for a solid species.
auto standardThermoModelSolid(real T, real P) -> StandardThermoProps
{
    StandardThermoProps props;
    props.G0  = 1.1 * (T*P)*(T*P);
    props.H0  = 1.2 * (T*P)*(T*P);
    props.V0  = 1.3 * (T*P)*(T*P);
    props.Cp0 = 1.4 * (T*P)*(T*P);
    props.Cv0 = 1.5 * (T*P)*(T*P);
    return props;
};

/// Return a mock Species object of aqueous type for test reasons.
auto createAqueousSpecies(String name)
{
    return Species(name)
        .withAggregateState(AggregateState::Aqueous)
        .withStandardThermoModel(standardThermoModelAqueous);
}

/// Return a mock Species object of gas type for test reasons.
auto createGaseousSpecies(String name)
{
    return Species(name)
        .withAggregateState(AggregateState::Gas)
        .withStandardThermoModel(standardThermoModelGaseous);
}

/// Return a mock Species object of solid type for test reasons.
auto createSolidSpecies(String name)
{
    return Species(name)
        .withAggregateState(AggregateState::Solid)
        .withStandardThermoModel(standardThermoModelSolid);
}

/// Return a mock ChemicalSystem object for test reasons.
auto createDatabase() -> Database
{
    Database db;

    // AQUEOUS SPECIES
    db.addSpecies( createAqueousSpecies("H2O(aq)"  ) );
    db.addSpecies( createAqueousSpecies("H+(aq)"   ) );
    db.addSpecies( createAqueousSpecies("OH-(aq)"  ) );
    db.addSpecies( createAqueousSpecies("H2(aq)"   ) );
    db.addSpecies( createAqueousSpecies("O2(aq)"   ) );
    db.addSpecies( createAqueousSpecies("Na+(aq)"  ) );
    db.addSpecies( createAqueousSpecies("Cl-(aq)"  ) );
    db.addSpecies( createAqueousSpecies("NaCl(aq)" ) );
    db.addSpecies( createAqueousSpecies("HCl(aq)"  ) );
    db.addSpecies( createAqueousSpecies("NaOH(aq)" ) );
    db.addSpecies( createAqueousSpecies("Ca++(aq)" ) );
    db.addSpecies( createAqueousSpecies("Mg++(aq)" ) );
    db.addSpecies( createAqueousSpecies("CO2(aq)"  ) );
    db.addSpecies( createAqueousSpecies("HCO3-(aq)") );
    db.addSpecies( createAqueousSpecies("CO3--(aq)") );
    db.addSpecies( createAqueousSpecies("CaCl2(aq)") );
    db.addSpecies( createAqueousSpecies("MgCl2(aq)") );
    db.addSpecies( createAqueousSpecies("SiO2(aq)" ) );

    // GASEOUS SPECIES
    db.addSpecies( createGaseousSpecies("CO2(g)" ) );
    db.addSpecies( createGaseousSpecies("O2(g)"  ) );
    db.addSpecies( createGaseousSpecies("H2(g)"  ) );
    db.addSpecies( createGaseousSpecies("H2O(g)" ) );
    db.addSpecies( createGaseousSpecies("CH4(g)" ) );
    db.addSpecies( createGaseousSpecies("CO(g)"  ) );

    // SOLID SPECIES
    db.addSpecies( createSolidSpecies("NaCl(s)"      ) );
    db.addSpecies( createSolidSpecies("CaCO3(s)"     ) );
    db.addSpecies( createSolidSpecies("MgCO3(s)"     ) );
    db.addSpecies( createSolidSpecies("CaMg(CO3)2(s)") );
    db.addSpecies( createSolidSpecies("SiO2(s)"      ) );

    return db;
}

} // namespace test

TEST_CASE("Testing Database class", "[Database]")
{
    Database db = test::createDatabase();

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::elements
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.elements().index("H") );
    REQUIRE_NOTHROW( db.elements().index("C") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::species
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.species().index("H2O(aq)") );
    REQUIRE_NOTHROW( db.species().index("CO2(g)")  );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::speciesWithAggregateState
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Aqueous).index("H2O(aq)") );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Aqueous).index("H+(aq)")  );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Gas).index("CO2(g)")      );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Solid).index("NaCl(s)")   );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::addSpecies
    //-------------------------------------------------------------------------
    db.addSpecies( Species("Fe(s)").withName("Iron")   );
    db.addSpecies( Species("Cu(s)").withName("Copper") );

    REQUIRE_NOTHROW( db.species().index("Iron")   );
    REQUIRE_NOTHROW( db.species().index("Copper") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::attachData
    //-------------------------------------------------------------------------
    db.attachData( String("SomeData") );

    REQUIRE( std::any_cast<String>(db.attachedData()) == "SomeData" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::clear
    //-------------------------------------------------------------------------
    db.clear();

    REQUIRE( db.elements().size() == 0 );
    REQUIRE( db.species().size()  == 0 );
}
