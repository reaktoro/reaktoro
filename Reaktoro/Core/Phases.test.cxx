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
    db.addSpecies( Species("H2O(aq)") );
    db.addSpecies( Species("H+") );
    db.addSpecies( Species("OH-") );
    db.addSpecies( Species("H2(aq)") );
    db.addSpecies( Species("O2(aq)") );
    db.addSpecies( Species("Na+") );
    db.addSpecies( Species("Cl-") );
    db.addSpecies( Species("NaCl(aq)") );
    db.addSpecies( Species("HCl(aq)") );
    db.addSpecies( Species("NaOH(aq)") );
    db.addSpecies( Species("Ca++") );
    db.addSpecies( Species("Mg++") );
    db.addSpecies( Species("CO2(aq)") );
    db.addSpecies( Species("HCO3-") );
    db.addSpecies( Species("CO3--") );
    db.addSpecies( Species("CaCl2(aq)") );
    db.addSpecies( Species("MgCl2(aq)") );
    db.addSpecies( Species("SiO2(aq)") );
    db.addSpecies( Species("C4H9OH(aq)").withTags("organic").withName("1-Butanol(aq)") );
    db.addSpecies( Species("C4H8(aq)").withTags("organic").withName("1-Butene(aq)") );
    db.addSpecies( Species("NaCl(s)").withName("Halite") );
    db.addSpecies( Species("CaCO3(s)").withName("Calcite").withTags("carbonate") );
    db.addSpecies( Species("MgCO3(s)").withName("Magnesite").withTags("carbonate") );
    db.addSpecies( Species("CaMg(CO3)2(s)").withName("Dolomite").withTags("carbonate") );
    db.addSpecies( Species("SiO2(s)").withName("Quartz") );
    db.addSpecies( Species("C(s)").withName("Graphite") );
    db.addSpecies( Species("CaO(s)").withName("Lime") );
    db.addSpecies( Species("N2(g)").withTags("inert") );
    db.addSpecies( Species("BaSO4(s)").withName("Barite").withTags("sulfate") );
    db.addSpecies( Species("SrSO4(s)").withName("Celestite").withTags("sulfate") );
    db.addSpecies( Species("PbSO4(s)").withName("Anglesite").withTags("sulfate") );
    db.addSpecies( Species("CaSO4(s)").withName("Anhydrite").withTags("sulfate") );
    db.addSpecies( Species("CO2(g)") );
    db.addSpecies( Species("O2(g)") );
    db.addSpecies( Species("H2(g)") );
    db.addSpecies( Species("H2O(g)") );
    db.addSpecies( Species("CH4(g)") );
    db.addSpecies( Species("CO(g)") );

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

        CHECK( phase.name() == "AqueousPhase" );
        CHECK( phase.stateOfMatter() == StateOfMatter::Liquid );
        CHECK( phase.species().size() == 3 );
        CHECK( phase.species(0).name() == "H2O(aq)" );
        CHECK( phase.species(1).name() == "H+"      );
        CHECK( phase.species(2).name() == "OH-"     );
        CHECK( phase.activityPropsFn() );
        CHECK( phase.idealActivityPropsFn() );
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

        CHECK( phase.name() == "AqueousPhase" );
        CHECK( phase.stateOfMatter() == StateOfMatter::Liquid );
        CHECK( phase.species().size() == 5 );
        CHECK( phase.species(0).name() == "H2O(aq)" );
        CHECK( phase.species(1).name() == "H+"      );
        CHECK( phase.species(2).name() == "OH-"     );
        CHECK( phase.species(3).name() == "H2(aq)"  );
        CHECK( phase.species(4).name() == "O2(aq)"  );
        CHECK( phase.activityPropsFn() );
        CHECK( phase.idealActivityPropsFn() );
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

        CHECK( phase.name() == "AqueousPhase" );
        CHECK( phase.stateOfMatter() == StateOfMatter::Liquid );
        CHECK( phase.species().size() == 20 );
        CHECK( phase.species(0).name()  == "H2O(aq)"       );
        CHECK( phase.species(1).name()  == "H+"            );
        CHECK( phase.species(2).name()  == "OH-"           );
        CHECK( phase.species(3).name()  == "H2(aq)"        );
        CHECK( phase.species(4).name()  == "O2(aq)"        );
        CHECK( phase.species(5).name()  == "Na+"           );
        CHECK( phase.species(6).name()  == "Cl-"           );
        CHECK( phase.species(7).name()  == "NaCl(aq)"      );
        CHECK( phase.species(8).name()  == "HCl(aq)"       );
        CHECK( phase.species(9).name()  == "NaOH(aq)"      );
        CHECK( phase.species(10).name() == "Ca++"          );
        CHECK( phase.species(11).name() == "Mg++"          );
        CHECK( phase.species(12).name() == "CO2(aq)"       );
        CHECK( phase.species(13).name() == "HCO3-"         );
        CHECK( phase.species(14).name() == "CO3--"         );
        CHECK( phase.species(15).name() == "CaCl2(aq)"     );
        CHECK( phase.species(16).name() == "MgCl2(aq)"     );
        CHECK( phase.species(17).name() == "SiO2(aq)"      );
        CHECK( phase.species(18).name() == "1-Butanol(aq)" );
        CHECK( phase.species(19).name() == "1-Butene(aq)"  );
        CHECK( phase.activityPropsFn() );
        CHECK( phase.idealActivityPropsFn() );
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

        CHECK( phase.name() == "AqueousPhase" );
        CHECK( phase.stateOfMatter() == StateOfMatter::Liquid );
        CHECK( phase.species().size() == 7 );
        CHECK( phase.species(0).name() == "H2O(aq)"  );
        CHECK( phase.species(1).name() == "H+"       );
        CHECK( phase.species(2).name() == "OH-"      );
        CHECK( phase.species(3).name() == "H2(aq)"   );
        CHECK( phase.species(4).name() == "O2(aq)"   );
        CHECK( phase.species(5).name() == "Cl-"      );
        CHECK( phase.species(6).name() == "HCl(aq)"  );
        CHECK( phase.activityPropsFn() );
        CHECK( phase.idealActivityPropsFn() );
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

        CHECK( phases.size() == 2 );

        CHECK( phases[0].name() == "Calcite" );
        CHECK( phases[0].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[0].aggregateState() == AggregateState::Solid );
        CHECK( phases[0].activityModel() );
        CHECK( phases[0].idealActivityModel() );

        CHECK( phases[1].name() == "Magnesite" );
        CHECK( phases[1].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[1].aggregateState() == AggregateState::Solid );
        CHECK( phases[1].activityModel() );
        CHECK( phases[1].idealActivityModel() );
    }

    SECTION("Testing GenericPhasesGenerator::GenericPhasesGenerator(const Speciate&)")
    {
        GenericPhasesGenerator generator(speciate("Ca Mg C O"));
        generator.setStateOfMatter(StateOfMatter::Solid);
        generator.setAggregateState(AggregateState::Solid);
        generator.setActivityModel(activitymodel);
        generator.setIdealActivityModel(activitymodel);

        Vec<GenericPhase> phases = generator.convert(db, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        CHECK( phases.size() == 5 );

        CHECK( phases[0].name() == "Calcite" );
        CHECK( phases[0].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[0].aggregateState() == AggregateState::Solid );
        CHECK( phases[0].activityModel() );
        CHECK( phases[0].idealActivityModel() );

        CHECK( phases[1].name() == "Magnesite" );
        CHECK( phases[1].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[1].aggregateState() == AggregateState::Solid );
        CHECK( phases[1].activityModel() );
        CHECK( phases[1].idealActivityModel() );

        CHECK( phases[2].name() == "Dolomite" );
        CHECK( phases[2].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[2].aggregateState() == AggregateState::Solid );
        CHECK( phases[2].activityModel() );
        CHECK( phases[2].idealActivityModel() );

        CHECK( phases[3].name() == "Graphite" );
        CHECK( phases[3].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[3].aggregateState() == AggregateState::Solid );
        CHECK( phases[3].activityModel() );
        CHECK( phases[3].idealActivityModel() );

        CHECK( phases[4].name() == "Lime" );
        CHECK( phases[4].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[4].aggregateState() == AggregateState::Solid );
        CHECK( phases[4].activityModel() );
        CHECK( phases[4].idealActivityModel() );
    }

    SECTION("Testing GenericPhasesGenerator::GenericPhasesGenerator() with all elements during phase conversion process")
    {
        GenericPhasesGenerator generator;
        generator.setStateOfMatter(StateOfMatter::Solid);
        generator.setAggregateState(AggregateState::Solid);
        generator.setActivityModel(activitymodel);
        generator.setIdealActivityModel(activitymodel);

        Vec<GenericPhase> phases = generator.convert(db, {"H", "O", "C", "Na", "Cl", "Ca", "Mg", "Si"});

        CHECK( phases.size() == 7 );

        CHECK( phases[0].name() == "Halite" );
        CHECK( phases[0].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[0].aggregateState() == AggregateState::Solid );
        CHECK( phases[0].activityModel() );
        CHECK( phases[0].idealActivityModel() );

        CHECK( phases[1].name() == "Calcite" );
        CHECK( phases[1].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[1].aggregateState() == AggregateState::Solid );
        CHECK( phases[1].activityModel() );
        CHECK( phases[1].idealActivityModel() );

        CHECK( phases[2].name() == "Magnesite" );
        CHECK( phases[2].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[2].aggregateState() == AggregateState::Solid );
        CHECK( phases[2].activityModel() );
        CHECK( phases[2].idealActivityModel() );

        CHECK( phases[3].name() == "Dolomite" );
        CHECK( phases[3].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[3].aggregateState() == AggregateState::Solid );
        CHECK( phases[3].activityModel() );
        CHECK( phases[3].idealActivityModel() );

        CHECK( phases[4].name() == "Quartz" );
        CHECK( phases[4].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[4].aggregateState() == AggregateState::Solid );
        CHECK( phases[4].activityModel() );
        CHECK( phases[4].idealActivityModel() );

        CHECK( phases[5].name() == "Graphite" );
        CHECK( phases[5].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[5].aggregateState() == AggregateState::Solid );
        CHECK( phases[5].activityModel() );
        CHECK( phases[5].idealActivityModel() );

        CHECK( phases[6].name() == "Lime" );
        CHECK( phases[6].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[6].aggregateState() == AggregateState::Solid );
        CHECK( phases[6].activityModel() );
        CHECK( phases[6].idealActivityModel() );
    }

    SECTION("Testing GenericPhasesGenerator::GenericPhasesGenerator() with Na and Cl elements during phase conversion process")
    {
        GenericPhasesGenerator generator;
        generator.setStateOfMatter(StateOfMatter::Solid);
        generator.setAggregateState(AggregateState::Solid);
        generator.setActivityModel(activitymodel);
        generator.setIdealActivityModel(activitymodel);

        Vec<GenericPhase> phases = generator.convert(db, {"Na", "Cl"});

        CHECK( phases.size() == 1 );

        CHECK( phases[0].name() == "Halite" );
        CHECK( phases[0].stateOfMatter() == StateOfMatter::Solid );
        CHECK( phases[0].aggregateState() == AggregateState::Solid );
        CHECK( phases[0].activityModel() );
        CHECK( phases[0].idealActivityModel() );
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: Phases
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================

    auto checkPhase = [&](Phase phase, String name, StringList species, StateOfMatter stateofmatter, AggregateState aggregatestate)
    {
        CHECK( phase.name() == name );
        CHECK( phase.stateOfMatter() == stateofmatter );
        CHECK( phase.aggregateState() == aggregatestate );
        CHECK( phase.activityPropsFn() );
        CHECK( phase.idealActivityPropsFn() );
        CHECK( phase.species().size() == species.size() );
        for(auto i = 0; i < species.size(); ++i)
        {
            INFO(str("phase name = ", name, ", species index = ", i));
            CHECK( phase.species(i).name() == species[i] );
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

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 1 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- Na+ Cl- HCO3- CO3-- CO2(aq)");
    }

    SECTION("Testing Phases with an aqueous solution using speciate together with a gaseous phase and a mineral")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("H O Na Cl C")) );
        phases.add( GaseousPhase("H2O(g) CO2(g)") );
        phases.add( MineralPhase("Halite") );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 3 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) CO2(aq) HCO3- CO3-- 1-Butanol(aq) 1-Butene(aq)");
        checkGaseousPhase(phasevec[1], "H2O(g) CO2(g)");
        checkMineralPhase(phasevec[2], "Halite");
    }

    SECTION("Testing Phases with an aqueous solution speciated automatically together with a gaseous phase and minerals")
    {
        Phases phases(db);

        phases.add( AqueousPhase() );
        phases.add( GaseousPhase(speciate("H O C")) );
        phases.add( MineralPhases("Halite Calcite Magnesite Dolomite Quartz") );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 7 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) Ca++ Mg++ CO2(aq) HCO3- CO3-- CaCl2(aq) MgCl2(aq) SiO2(aq) 1-Butanol(aq) 1-Butene(aq)");
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

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 6 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) Ca++ CO2(aq) HCO3- CO3-- CaCl2(aq) 1-Butanol(aq) 1-Butene(aq)");
        checkGaseousPhase(phasevec[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)");
        checkMineralPhase(phasevec[2], "Halite");
        checkMineralPhase(phasevec[3], "Calcite");
        checkMineralPhase(phasevec[4], "Graphite");
        checkMineralPhase(phasevec[5], "Lime");
    }

    SECTION("Testing Phases with minerals collected automatically from given aqueous species names")
    {
        Phases phases(db);

        phases.add( AqueousPhase("H2O(aq) H+ OH- Na+ Cl- CO2(aq) Ca++") );
        phases.add( MineralPhases() );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 5 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- Na+ Cl- CO2(aq) Ca++");
        checkMineralPhase(phasevec[1], "Halite");
        checkMineralPhase(phasevec[2], "Calcite");
        checkMineralPhase(phasevec[3], "Graphite");
        checkMineralPhase(phasevec[4], "Lime");
    }

    SECTION("Testing Phases with no minerals collected because none exist with H, O, C elements")
    {
        Phases phases(db);

        phases.add( AqueousPhase() );
        phases.add( GaseousPhase(speciate("H O C")) );
        phases.add( MineralPhases() ); // no mineral should be collected because there is none with elements {H, O, C}

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 3 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) CO2(aq) HCO3- CO3-- 1-Butanol(aq) 1-Butene(aq)");
        checkGaseousPhase(phasevec[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)");
        checkMineralPhase(phasevec[2], "Graphite");
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: AqueousPhase with provided speciates and tags, so that species possessing them are excluded
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing AqueousPhase::AqueousPhase(Speciate, Exclude)")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("H O C Na Cl")) );
        phases.add( GaseousPhase("CO2(g)") );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 2 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) CO2(aq) HCO3- CO3-- 1-Butanol(aq) 1-Butene(aq)");
        checkGaseousPhase(phasevec[1], "CO2(g)");
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: AqueousPhase with provided speciates and tags, so that species possessing them are excluded
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing AqueousPhase::AqueousPhase(Speciate, Exclude)")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("H O C Na Cl"), exclude("organic")) );
        phases.add( GaseousPhase("CO2(g)") );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 2 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) CO2(aq) HCO3- CO3--");
        checkGaseousPhase(phasevec[1], "CO2(g)");
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: AqueousPhase with provided tags, so that species possessing them are excluded
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing AqueousPhase::AqueousPhase(Exclude)")
    {
        Phases phases(db);

        phases.add( AqueousPhase(exclude("organic")) );
        phases.add( GaseousPhase("CO2(g)") );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 2 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq)");
        checkGaseousPhase(phasevec[1], "CO2(g)");
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: MineralPhases with provided tags, so that species possessing them are excluded
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing MineralPhases::MineralPhases(Exclude)")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("H O C"), exclude("organic")) );
        phases.add( GaseousPhase("H2O(g) CO2(g)") );
        phases.add( MineralPhases(exclude("carbonate")) );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 3 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) CO2(aq) HCO3- CO3--");
        checkGaseousPhase(phasevec[1], "H2O(g) CO2(g)");
        checkMineralPhase(phasevec[2], "Graphite");
    }
//
    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: MineralPhases with provided speciate symbols and tags, so that species possessing them are excluded
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing MineralPhases::MineralPhases(Speciate, Exclude)")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("H O"), exclude("organic")) );
        phases.add( MineralPhases(speciate("C Ca O"), exclude("carbonate")) );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 3 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq)");
        checkMineralPhase(phasevec[1], "Graphite");
        checkMineralPhase(phasevec[2], "Lime");
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: GaseousPhases with provided speciates and tags, so that species possessing them are excluded
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing GaseousPhases::GaseousPhases(Speciate, Exclude)")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("H O C"), exclude("organic")) );
        phases.add( GaseousPhase(speciate("H O C"), exclude("inert")) );
        phases.add( MineralPhases(exclude("carbonate")) );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 3 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) CO2(aq) HCO3- CO3--");
        checkGaseousPhase(phasevec[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)");
        checkMineralPhase(phasevec[2], "Graphite");
    }

    //=================================================================================================================
    //-----------------------------------------------------------------------------------------------------------------
    // TESTING CLASS: GaseousPhases with provided tags, so that species possessing them are excluded
    //-----------------------------------------------------------------------------------------------------------------
    //=================================================================================================================
    SECTION("Testing GaseousPhases::GaseousPhases(Exclude)")
    {
        Phases phases(db);

        phases.add( AqueousPhase(speciate("H O C Na Cl"), exclude("organic")) );
        phases.add( GaseousPhase(exclude("inert")) );
        phases.add( MineralPhases(exclude("carbonate")) );

        Vec<Phase> phasevec = phases.convert();

        CHECK( phasevec.size() == 4 );

        checkAqueousPhase(phasevec[0], "H2O(aq) H+ OH- H2(aq) O2(aq) Na+ Cl- NaCl(aq) HCl(aq) NaOH(aq) CO2(aq) HCO3- CO3--");
        checkGaseousPhase(phasevec[1], "CO2(g) O2(g) H2(g) H2O(g) CH4(g) CO(g)");
        checkMineralPhase(phasevec[2], "Halite");
        checkMineralPhase(phasevec[3], "Graphite");
    }
}
