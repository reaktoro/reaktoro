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
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ThermoPropsPhase.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Phases", "[Phases]")
{
    StandardThermoPropsFn standard_thermo_props_fn = [](real T, real P, const Species& species)
    {
        return StandardThermoProps{};
    };

    ActivityModel activity_model = [](const SpeciesList& species)
    {
        ActivityPropsFn fn = [](ActivityProps props, real T, real P, ArrayXrConstRef x) {};
        return fn;
    };

    Database db;
    db.addSpecies( Species("H2O(aq)")                            );
    db.addSpecies( Species("H+")                                 );
    db.addSpecies( Species("OH-")                                );
    db.addSpecies( Species("H2(aq)")                             );
    db.addSpecies( Species("O2(aq)")                             );
    db.addSpecies( Species("Na+")                                );
    db.addSpecies( Species("Cl-")                                );
    db.addSpecies( Species("NaCl(aq)")                           );
    db.addSpecies( Species("HCl(aq)")                            );
    db.addSpecies( Species("NaOH(aq)")                           );
    db.addSpecies( Species("Ca++")                               );
    db.addSpecies( Species("Mg++")                               );
    db.addSpecies( Species("CaCl2(aq)")                          );
    db.addSpecies( Species("MgCl2(aq)")                          );
    db.addSpecies( Species("SiO2(aq)")                           );
    db.addSpecies( Species("NaCl(s)").withName("Halite")         );
    db.addSpecies( Species("CaCO3(s)").withName("Calcite")       );
    db.addSpecies( Species("MgCO3(s)").withName("Magnesite")     );
    db.addSpecies( Species("CaMg(CO3)2(s)").withName("Dolomite") );
    db.addSpecies( Species("SiO2(s)").withName("Quartz")         );

    ThermoEngine engine(db, standard_thermo_props_fn);

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: GenericPhase
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing GenericPhase::GenericPhase(const StringList&)")
    {
        GenericPhase genericphase("H2O(aq) H+ OH-");
        genericphase.setName("AqueousSolution");
        genericphase.setStateOfMatter(StateOfMatter::Liquid);
        genericphase.setAggregateState(AggregateState::Aqueous);
        genericphase.setActivityModel(activity_model);

        Phase phase = genericphase.convert(engine, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phase.name() == "AqueousSolution" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 3 );
        REQUIRE( phase.species(0).name() == "H2O(aq)" );
        REQUIRE( phase.species(1).name() == "H+"      );
        REQUIRE( phase.species(2).name() == "OH-"     );
        REQUIRE( phase.standardThermoPropsFn() );
        REQUIRE( phase.activityPropsFn() );
    }

    SECTION("Testing GenericPhase::GenericPhase(const Speciate&)")
    {
        GenericPhase genericphase(speciate("H O"));
        genericphase.setName("AqueousSolution");
        genericphase.setStateOfMatter(StateOfMatter::Liquid);
        genericphase.setAggregateState(AggregateState::Aqueous);
        genericphase.setActivityModel(activity_model);

        Phase phase = genericphase.convert(engine, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phase.name() == "AqueousSolution" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 5 );
        REQUIRE( phase.species(0).name() == "H2O(aq)" );
        REQUIRE( phase.species(1).name() == "H+"      );
        REQUIRE( phase.species(2).name() == "OH-"     );
        REQUIRE( phase.species(3).name() == "H2(aq)"  );
        REQUIRE( phase.species(4).name() == "O2(aq)"  );
        REQUIRE( phase.standardThermoPropsFn() );
        REQUIRE( phase.activityPropsFn() );
    }

    SECTION("Testing GenericPhase::GenericPhase() with all elements during phase conversion process")
    {
        GenericPhase genericphase;
        genericphase.setName("AqueousSolution");
        genericphase.setStateOfMatter(StateOfMatter::Liquid);
        genericphase.setAggregateState(AggregateState::Aqueous);
        genericphase.setActivityModel(activity_model);

        Phase phase = genericphase.convert(engine, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phase.name() == "AqueousSolution" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 15 );
        REQUIRE( phase.species(0).name()  == "H2O(aq)"   );
        REQUIRE( phase.species(1).name()  == "H+"        );
        REQUIRE( phase.species(2).name()  == "OH-"       );
        REQUIRE( phase.species(3).name()  == "H2(aq)"    );
        REQUIRE( phase.species(4).name()  == "O2(aq)"    );
        REQUIRE( phase.species(5).name()  == "Na+"       );
        REQUIRE( phase.species(6).name()  == "Cl-"       );
        REQUIRE( phase.species(7).name()  == "NaCl(aq)"  );
        REQUIRE( phase.species(8).name()  == "HCl(aq)"   );
        REQUIRE( phase.species(9).name()  == "NaOH(aq)"  );
        REQUIRE( phase.species(10).name() == "Ca++"      );
        REQUIRE( phase.species(11).name() == "Mg++"      );
        REQUIRE( phase.species(12).name() == "CaCl2(aq)" );
        REQUIRE( phase.species(13).name() == "MgCl2(aq)" );
        REQUIRE( phase.species(14).name() == "SiO2(aq)"  );
        REQUIRE( phase.standardThermoPropsFn() );
        REQUIRE( phase.activityPropsFn() );
    }

    SECTION("Testing GenericPhase::GenericPhase() with H, O, and Cl elements during phase conversion process")
    {
        GenericPhase genericphase;
        genericphase.setName("AqueousSolution");
        genericphase.setStateOfMatter(StateOfMatter::Liquid);
        genericphase.setAggregateState(AggregateState::Aqueous);
        genericphase.setActivityModel(activity_model);

        Phase phase = genericphase.convert(engine, {"H", "O", "Cl"});

        REQUIRE( phase.name() == "AqueousSolution" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 7 );
        REQUIRE( phase.species(0).name() == "H2O(aq)"  );
        REQUIRE( phase.species(1).name() == "H+"       );
        REQUIRE( phase.species(2).name() == "OH-"      );
        REQUIRE( phase.species(3).name() == "H2(aq)"   );
        REQUIRE( phase.species(4).name() == "O2(aq)"   );
        REQUIRE( phase.species(5).name() == "Cl-"      );
        REQUIRE( phase.species(6).name() == "HCl(aq)"  );
        REQUIRE( phase.standardThermoPropsFn() );
        REQUIRE( phase.activityPropsFn() );
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: GenericPhases
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing GenericPhases::GenericPhases(const StringList&)")
    {
        GenericPhases genericphases("Calcite Magnesite");
        genericphases.setStateOfMatter(StateOfMatter::Solid);
        genericphases.setAggregateState(AggregateState::Solid);
        genericphases.setActivityModel(activity_model);

        Vec<GenericPhase> phases = genericphases.convert(engine, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phases.size() == 2 );

        REQUIRE( phases[0].name() == "Calcite" );
        REQUIRE( phases[0].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[0].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[0].activityModel() );

        REQUIRE( phases[1].name() == "Magnesite" );
        REQUIRE( phases[1].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[1].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[1].activityModel() );
    }

    SECTION("Testing GenericPhases::GenericPhases(const Speciate&)")
    {
        GenericPhases genericphases(speciate("Ca Mg C O"));
        genericphases.setStateOfMatter(StateOfMatter::Solid);
        genericphases.setAggregateState(AggregateState::Solid);
        genericphases.setActivityModel(activity_model);

        Vec<GenericPhase> phases = genericphases.convert(engine, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phases.size() == 3 );

        REQUIRE( phases[0].name() == "Calcite" );
        REQUIRE( phases[0].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[0].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[0].activityModel() );

        REQUIRE( phases[1].name() == "Magnesite" );
        REQUIRE( phases[1].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[1].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[1].activityModel() );

        REQUIRE( phases[2].name() == "Dolomite" );
        REQUIRE( phases[2].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[2].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[2].activityModel() );
    }

    SECTION("Testing GenericPhases::GenericPhases() with all elements during phase conversion process")
    {
        GenericPhases genericphases;
        genericphases.setStateOfMatter(StateOfMatter::Solid);
        genericphases.setAggregateState(AggregateState::Solid);
        genericphases.setActivityModel(activity_model);

        Vec<GenericPhase> phases = genericphases.convert(engine, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phases.size() == 5 );

        REQUIRE( phases[0].name() == "Halite" );
        REQUIRE( phases[0].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[0].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[0].activityModel() );

        REQUIRE( phases[1].name() == "Calcite" );
        REQUIRE( phases[1].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[1].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[1].activityModel() );

        REQUIRE( phases[2].name() == "Magnesite" );
        REQUIRE( phases[2].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[2].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[2].activityModel() );

        REQUIRE( phases[3].name() == "Dolomite" );
        REQUIRE( phases[3].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[3].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[3].activityModel() );

        REQUIRE( phases[4].name() == "Quartz" );
        REQUIRE( phases[4].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[4].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[4].activityModel() );
    }

    SECTION("Testing GenericPhases::GenericPhases() with Na and Cl elements during phase conversion process")
    {
        GenericPhases genericphases;
        genericphases.setStateOfMatter(StateOfMatter::Solid);
        genericphases.setAggregateState(AggregateState::Solid);
        genericphases.setActivityModel(activity_model);

        Vec<GenericPhase> phases = genericphases.convert(engine, {"Na", "Cl"});

        REQUIRE( phases.size() == 1 );

        REQUIRE( phases[0].name() == "Halite" );
        REQUIRE( phases[0].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[0].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[0].activityModel() );
    }
}
