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
using namespace Reaktoro;

TEST_CASE("Testing Phases", "[Phases]")
{
    ActivityModel activitymodel = [](const SpeciesList& species)
    {
        ActivityPropsFn fn = [](ActivityPropsRef props, ActivityArgs args) {};
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
    db.addSpecies( Species("CO2(aq)")                            );
    db.addSpecies( Species("HCO3-")                              );
    db.addSpecies( Species("CO3--")                              );
    db.addSpecies( Species("CaCl2(aq)")                          );
    db.addSpecies( Species("MgCl2(aq)")                          );
    db.addSpecies( Species("SiO2(aq)")                           );
    db.addSpecies( Species("NaCl(s)").withName("Halite")         );
    db.addSpecies( Species("CaCO3(s)").withName("Calcite")       );
    db.addSpecies( Species("MgCO3(s)").withName("Magnesite")     );
    db.addSpecies( Species("CaMg(CO3)2(s)").withName("Dolomite") );
    db.addSpecies( Species("SiO2(s)").withName("Quartz")         );
    db.addSpecies( Species("CO2(g)")                             );
    db.addSpecies( Species("O2(g)")                              );
    db.addSpecies( Species("H2(g)")                              );
    db.addSpecies( Species("H2O(g)")                             );
    db.addSpecies( Species("CH4(g)")                             );
    db.addSpecies( Species("CO(g)")                              );

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: GenericPhase
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing GenericPhase::GenericPhase(const StringList&)")
    {
        GenericPhase genericphase("H2O(aq) H+ OH-");
        genericphase.setName("AqueousPhase");
        genericphase.setStateOfMatter(StateOfMatter::Liquid);
        genericphase.setAggregateState(AggregateState::Aqueous);
        genericphase.setActivityModel(activitymodel);
        genericphase.setIdealActivityModel(activitymodel);

        Phase phase = genericphase.convert(db, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phase.name() == "AqueousPhase" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 3 );
        REQUIRE( phase.species(0).name() == "H2O(aq)" );
        REQUIRE( phase.species(1).name() == "H+"      );
        REQUIRE( phase.species(2).name() == "OH-"     );
        REQUIRE( phase.activityPropsFn() );
        REQUIRE( phase.idealActivityPropsFn() );
    }

    SECTION("Testing GenericPhase::GenericPhase(const Speciate&)")
    {
        GenericPhase genericphase(speciate("H O"));
        genericphase.setName("AqueousPhase");
        genericphase.setStateOfMatter(StateOfMatter::Liquid);
        genericphase.setAggregateState(AggregateState::Aqueous);
        genericphase.setActivityModel(activitymodel);
        genericphase.setIdealActivityModel(activitymodel);

        Phase phase = genericphase.convert(db, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phase.name() == "AqueousPhase" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 5 );
        REQUIRE( phase.species(0).name() == "H2O(aq)" );
        REQUIRE( phase.species(1).name() == "H+"      );
        REQUIRE( phase.species(2).name() == "OH-"     );
        REQUIRE( phase.species(3).name() == "H2(aq)"  );
        REQUIRE( phase.species(4).name() == "O2(aq)"  );
        REQUIRE( phase.activityPropsFn() );
        REQUIRE( phase.idealActivityPropsFn() );
    }

    SECTION("Testing GenericPhase::GenericPhase() with all elements during phase conversion process")
    {
        GenericPhase genericphase;
        genericphase.setName("AqueousPhase");
        genericphase.setStateOfMatter(StateOfMatter::Liquid);
        genericphase.setAggregateState(AggregateState::Aqueous);
        genericphase.setActivityModel(activitymodel);
        genericphase.setIdealActivityModel(activitymodel);

        Phase phase = genericphase.convert(db, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phase.name() == "AqueousPhase" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 18 );
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
        REQUIRE( phase.species(12).name() == "CO2(aq)"   );
        REQUIRE( phase.species(13).name() == "HCO3-"     );
        REQUIRE( phase.species(14).name() == "CO3--"     );
        REQUIRE( phase.species(15).name() == "CaCl2(aq)" );
        REQUIRE( phase.species(16).name() == "MgCl2(aq)" );
        REQUIRE( phase.species(17).name() == "SiO2(aq)"  );
        REQUIRE( phase.activityPropsFn() );
        REQUIRE( phase.idealActivityPropsFn() );
    }

    SECTION("Testing GenericPhase::GenericPhase() with H, O, and Cl elements during phase conversion process")
    {
        GenericPhase genericphase;
        genericphase.setName("AqueousPhase");
        genericphase.setStateOfMatter(StateOfMatter::Liquid);
        genericphase.setAggregateState(AggregateState::Aqueous);
        genericphase.setActivityModel(activitymodel);
        genericphase.setIdealActivityModel(activitymodel);

        Phase phase = genericphase.convert(db, {"H", "O", "Cl"});

        REQUIRE( phase.name() == "AqueousPhase" );
        REQUIRE( phase.stateOfMatter() == StateOfMatter::Liquid );
        REQUIRE( phase.species().size() == 7 );
        REQUIRE( phase.species(0).name() == "H2O(aq)"  );
        REQUIRE( phase.species(1).name() == "H+"       );
        REQUIRE( phase.species(2).name() == "OH-"      );
        REQUIRE( phase.species(3).name() == "H2(aq)"   );
        REQUIRE( phase.species(4).name() == "O2(aq)"   );
        REQUIRE( phase.species(5).name() == "Cl-"      );
        REQUIRE( phase.species(6).name() == "HCl(aq)"  );
        REQUIRE( phase.activityPropsFn() );
        REQUIRE( phase.idealActivityPropsFn() );
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: GenericPhasesGenerator
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing GenericPhasesGenerator::GenericPhasesGenerator(const StringList&)")
    {
        GenericPhasesGenerator generator("Calcite Magnesite");
        generator.setStateOfMatter(StateOfMatter::Solid);
        generator.setAggregateState(AggregateState::Solid);
        generator.setActivityModel(activitymodel);
        generator.setIdealActivityModel(activitymodel);

        Vec<GenericPhase> phases = generator.convert(db, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phases.size() == 2 );

        REQUIRE( phases[0].name() == "Calcite" );
        REQUIRE( phases[0].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[0].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[0].activityModel() );
        REQUIRE( phases[0].idealActivityModel() );

        REQUIRE( phases[1].name() == "Magnesite" );
        REQUIRE( phases[1].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[1].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[1].activityModel() );
        REQUIRE( phases[1].idealActivityModel() );
    }

    SECTION("Testing GenericPhasesGenerator::GenericPhasesGenerator(const Speciate&)")
    {
        GenericPhasesGenerator generator(speciate("Ca Mg C O"));
        generator.setStateOfMatter(StateOfMatter::Solid);
        generator.setAggregateState(AggregateState::Solid);
        generator.setActivityModel(activitymodel);
        generator.setIdealActivityModel(activitymodel);

        Vec<GenericPhase> phases = generator.convert(db, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phases.size() == 3 );

        REQUIRE( phases[0].name() == "Calcite" );
        REQUIRE( phases[0].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[0].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[0].activityModel() );
        REQUIRE( phases[0].idealActivityModel() );

        REQUIRE( phases[1].name() == "Magnesite" );
        REQUIRE( phases[1].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[1].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[1].activityModel() );
        REQUIRE( phases[1].idealActivityModel() );

        REQUIRE( phases[2].name() == "Dolomite" );
        REQUIRE( phases[2].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[2].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[2].activityModel() );
        REQUIRE( phases[2].idealActivityModel() );
    }

    SECTION("Testing GenericPhasesGenerator::GenericPhasesGenerator() with all elements during phase conversion process")
    {
        GenericPhasesGenerator generator;
        generator.setStateOfMatter(StateOfMatter::Solid);
        generator.setAggregateState(AggregateState::Solid);
        generator.setActivityModel(activitymodel);
        generator.setIdealActivityModel(activitymodel);

        Vec<GenericPhase> phases = generator.convert(db, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        REQUIRE( phases.size() == 5 );

        REQUIRE( phases[0].name() == "Halite" );
        REQUIRE( phases[0].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[0].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[0].activityModel() );
        REQUIRE( phases[0].idealActivityModel() );

        REQUIRE( phases[1].name() == "Calcite" );
        REQUIRE( phases[1].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[1].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[1].activityModel() );
        REQUIRE( phases[1].idealActivityModel() );

        REQUIRE( phases[2].name() == "Magnesite" );
        REQUIRE( phases[2].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[2].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[2].activityModel() );
        REQUIRE( phases[2].idealActivityModel() );

        REQUIRE( phases[3].name() == "Dolomite" );
        REQUIRE( phases[3].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[3].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[3].activityModel() );
        REQUIRE( phases[3].idealActivityModel() );

        REQUIRE( phases[4].name() == "Quartz" );
        REQUIRE( phases[4].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[4].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[4].activityModel() );
        REQUIRE( phases[4].idealActivityModel() );
    }

    SECTION("Testing GenericPhasesGenerator::GenericPhasesGenerator() with Na and Cl elements during phase conversion process")
    {
        GenericPhasesGenerator generator;
        generator.setStateOfMatter(StateOfMatter::Solid);
        generator.setAggregateState(AggregateState::Solid);
        generator.setActivityModel(activitymodel);
        generator.setIdealActivityModel(activitymodel);

        Vec<GenericPhase> phases = generator.convert(db, {"Na", "Cl"});

        REQUIRE( phases.size() == 1 );

        REQUIRE( phases[0].name() == "Halite" );
        REQUIRE( phases[0].stateOfMatter() == StateOfMatter::Solid );
        REQUIRE( phases[0].aggregateState() == AggregateState::Solid );
        REQUIRE( phases[0].activityModel() );
        REQUIRE( phases[0].idealActivityModel() );
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: Phases
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================

    auto checkPhase = [&](Phase phase, String name, StringList species, StateOfMatter stateofmatter, AggregateState aggregatestate)
    {
        REQUIRE( phase.name() == name );
        REQUIRE( phase.stateOfMatter() == stateofmatter );
        REQUIRE( phase.aggregateState() == aggregatestate );
        REQUIRE( phase.activityPropsFn() );
        REQUIRE( phase.idealActivityPropsFn() );
        REQUIRE( phase.species().size() == species.size() );
        for(auto i = 0; i < species.size(); ++i)
        {
            INFO(str("phase name = ", name, ", species index = ", i));
            REQUIRE( phase.species(i).name() == species[i] );
        }
    };

    auto checkAqueousPhase = [&](Phase phase, StringList species)
    {
        checkPhase(phase, "AqueousPhase", species, StateOfMatter::Liquid, AggregateState::Aqueous);
    };

    auto checkGaseousPhase = [&](Phase phase, StringList species)
    {
        checkPhase(phase, "GaseousPhase", species, StateOfMatter::Gas, AggregateState::Gas);
    };

    auto checkMineralPhase = [&](Phase phase, StringList species)
    {
        checkPhase(phase, species[0], species, StateOfMatter::Solid, AggregateState::Solid);
    };

    SECTION("Testing Phases with an aqueous solution only")
    {
        Phases phases(db);
        phases.add( AqueousPhase("H2O(aq) H+ OH- Na+ Cl- HCO3- CO3-- CO2(aq)") );

        Vec<Phase> phasevec = phases;

        REQUIRE( phasevec.size() == 1 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- Na+ Cl- HCO3- CO3-- CO2(aq)");
    }

    SECTION("Testing Phases with an aqueous solution using speciate together with a gaseous phase and a mineral")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("H O Na Cl C")) );
        phases.add( GaseousPhase("H2O(g) CO2(g)") );
        phases.add( MineralPhase("Halite") );

        Vec<Phase> phasevec = phases;

        REQUIRE( phasevec.size() == 3 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) CO2(aq) HCO3- CO3--");
        checkGaseousPhase(phasevec[1], "H2O(g) CO2(g)");
        checkMineralPhase(phasevec[2], "Halite");
    }

    SECTION("Testing Phases with an aqueous solution speciated automatically together with a gaseous phase and minerals")
    {
        Phases phases(db);

        phases.add( AqueousPhase() );
        phases.add( GaseousPhase(speciate("H O C")) );
        phases.add( MineralPhases("Halite Calcite Magnesite Dolomite Quartz") );

        Vec<Phase> phasevec = phases;

        REQUIRE( phasevec.size() == 7 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) Ca++ Mg++ CO2(aq) HCO3- CO3-- CaCl2(aq) MgCl2(aq) SiO2(aq)");
        checkGaseousPhase(phasevec[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)");
        checkMineralPhase(phasevec[2], "Halite");
        checkMineralPhase(phasevec[3], "Calcite");
        checkMineralPhase(phasevec[4], "Magnesite");
        checkMineralPhase(phasevec[5], "Dolomite");
        checkMineralPhase(phasevec[6], "Quartz");
    }

    SECTION("Testing Phases with minerals collected automatically together with aqueous and gaseous phases")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("Na Cl C Ca")) );
        phases.add( GaseousPhase(speciate("H O C")) );
        phases.add( MineralPhases() );

        Vec<Phase> phasevec = phases;

        REQUIRE( phasevec.size() == 4 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) Ca++ CO2(aq) HCO3- CO3-- CaCl2(aq)");
        checkGaseousPhase(phasevec[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)");
        checkMineralPhase(phasevec[2], "Halite");
        checkMineralPhase(phasevec[3], "Calcite");
    }

    SECTION("Testing Phases with minerals collected automatically from given aqueous species names")
    {
        Phases phases(db);

        phases.add( AqueousPhase("H2O(aq) H+ OH- Na+ Cl- CO2(aq) Ca++") );
        phases.add( MineralPhases() );

        Vec<Phase> phasevec = phases;

        REQUIRE( phasevec.size() == 3 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- Na+ Cl- CO2(aq) Ca++");
        checkMineralPhase(phasevec[1], "Halite");
        checkMineralPhase(phasevec[2], "Calcite");
    }

    SECTION("Testing Phases with no minerals collected because none exist with H, O, C elements")
    {
        Phases phases(db);

        phases.add( AqueousPhase() );
        phases.add( GaseousPhase(speciate("H O C")) );
        phases.add( MineralPhases() ); // no mineral should be collected because there is none with elements {H, O, C}

        Vec<Phase> phasevec = phases;

        REQUIRE( phasevec.size() == 2 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) CO2(aq) HCO3- CO3--");
        checkGaseousPhase(phasevec[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)");
    }
}
