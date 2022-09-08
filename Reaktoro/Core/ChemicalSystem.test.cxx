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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Utils.hpp>
using namespace Reaktoro;

namespace test {

/// Return a mock Database object for test reasons (imported).
extern auto createDatabase() -> Database;

/// Return mock activity properties for an aqueous solution.
auto activityModelAqueous(ActivityPropsRef props, ActivityModelArgs args)
{
    const auto [T, P, x] = args;
    props.Vx  = 0.1 * log(P/T);
    props.VxT = 0.2 * log(P/T);
    props.VxP = 0.3 * log(P/T);
    props.Gx  = 0.4 * log(P/T);
    props.Hx  = 0.5 * log(P/T);
    props.Cpx = 0.6 * log(P/T);
    props.ln_g = 0.8 * x;
    props.ln_a = 0.9 * x;
};

/// Return mock activity properties for a gaseous solution.
auto activityModelGaseous(ActivityPropsRef props, ActivityModelArgs args)
{
    const auto [T, P, x] = args;
    props.Vx  = 1.0 * log(P/T);
    props.VxT = 2.0 * log(P/T);
    props.VxP = 3.0 * log(P/T);
    props.Gx  = 4.0 * log(P/T);
    props.Hx  = 5.0 * log(P/T);
    props.Cpx = 6.0 * log(P/T);
    props.ln_g = 8.0 * x;
    props.ln_a = 9.0 * x;
};

/// Return mock activity properties for a solid solution.
auto activityModelSolid(ActivityPropsRef props, ActivityModelArgs args)
{
    const auto [T, P, x] = args;
    props.Vx  = 1.1 * log(P/T);
    props.VxT = 2.1 * log(P/T);
    props.VxP = 3.1 * log(P/T);
    props.Gx  = 4.1 * log(P/T);
    props.Hx  = 5.1 * log(P/T);
    props.Cpx = 6.1 * log(P/T);
    props.ln_g = 8.1 * x;
    props.ln_a = 9.1 * x;
};

/// Return a Phase object representing an aqueous phase using mock thermodynamic models.
auto createAqueousPhase(const Database& db) -> Phase
{
    Phase phase;
    phase = phase.withName("AqueousPhase");
    phase = phase.withSpecies(db.speciesWithAggregateState(AggregateState::Aqueous));
    phase = phase.withStateOfMatter(StateOfMatter::Liquid);
    phase = phase.withActivityModel(activityModelAqueous);
    phase = phase.withIdealActivityModel(ActivityModelIdealAqueous()(phase.species()));
    return phase;
}

/// Return a Phase object representing a gaseous phase using mock thermodynamic models.
auto createGaseousPhase(const Database& db) -> Phase
{
    Phase phase;
    phase = phase.withName("GaseousPhase");
    phase = phase.withSpecies(db.speciesWithAggregateState(AggregateState::Gas));
    phase = phase.withStateOfMatter(StateOfMatter::Gas);
    phase = phase.withActivityModel(activityModelGaseous);
    phase = phase.withIdealActivityModel(ActivityModelIdealGas()(phase.species()));
    return phase;
}

/// Return a Phase object representing a mineral phase using mock thermodynamic models.
auto createMineralPhase(const Database& db, String phasename, String mineral) -> Phase
{
    Phase phase;
    phase = phase.withName(phasename);
    phase = phase.withSpecies({ db.species().get(mineral) });
    phase = phase.withStateOfMatter(StateOfMatter::Solid);
    phase = phase.withActivityModel(activityModelSolid);
    phase = phase.withIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)(phase.species()));
    return phase;
}

/// Return a mock ChemicalSystem object for test reasons.
auto createChemicalSystem() -> ChemicalSystem
{
    // Create the Database object for the ChemicalSystem
    Database db = test::createDatabase();

    // Create the Phase objects for the ChemicalSystem
    const PhaseList phases = {
        createAqueousPhase(db),
        createGaseousPhase(db),
        createMineralPhase(db, "Halite", "NaCl(s)"),
        createMineralPhase(db, "Calcite", "CaCO3(s)"),
        createMineralPhase(db, "Magnesite", "MgCO3(s)"),
        createMineralPhase(db, "Quartz", "SiO2(s)"),
    };

    const ReactionList reactions = {
        db.reaction("H2O(aq) = H+(aq) + OH-(aq)"),
        db.reaction("H2O(aq) = H2(aq) + 0.5*O2(aq)"),
        db.reaction("O2(g) = O2(aq)"),
        db.reaction("NaCl(s) = Na+(aq) + Cl-(aq)"),
        db.reaction("CaCO3(s)"),
    };

    return ChemicalSystem(db, phases, reactions);
}

} // namespace test

TEST_CASE("Testing ChemicalSystem class", "[ChemicalSystem]")
{
    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: ChemicalSystem::ChemicalSystem(phases)
    //-------------------------------------------------------------------------
    ChemicalSystem system = test::createChemicalSystem();

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::id()
    //-------------------------------------------------------------------------
    CHECK( ChemicalSystem().id() == system.id() + 1 ); // new ChemicalSystem object has id = id0 + 1 and original system continues to have id = id0
    CHECK( ChemicalSystem().id() == system.id() + 2 ); // new ChemicalSystem object has id = id0 + 2 and original system continues to have id = id0

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::database()
    //-------------------------------------------------------------------------
    CHECK( system.database().species().size() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::phases()
    //-------------------------------------------------------------------------
    CHECK( system.phases().size() == 6 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::phase(i)
    //-------------------------------------------------------------------------
    CHECK( system.phase(0).name() == "AqueousPhase" );
    CHECK( system.phase(1).name() == "GaseousPhase" );
    CHECK( system.phase(2).name() == "Halite"       );
    CHECK( system.phase(3).name() == "Calcite"      );
    CHECK( system.phase(4).name() == "Magnesite"    );
    CHECK( system.phase(5).name() == "Quartz"       );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::species()
    //-------------------------------------------------------------------------
    CHECK( system.species().size() == 29 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::species(i)
    //-------------------------------------------------------------------------
    CHECK( system.species(0).name()  == "H2O(aq)"   );
    CHECK( system.species(1).name()  == "H+(aq)"    );
    CHECK( system.species(2).name()  == "OH-(aq)"   );
    CHECK( system.species(3).name()  == "H2(aq)"    );
    CHECK( system.species(4).name()  == "O2(aq)"    );
    CHECK( system.species(5).name()  == "Na+(aq)"   );
    CHECK( system.species(6).name()  == "Cl-(aq)"   );
    CHECK( system.species(7).name()  == "NaCl(aq)"  );
    CHECK( system.species(8).name()  == "HCl(aq)"   );
    CHECK( system.species(9).name()  == "NaOH(aq)"  );
    CHECK( system.species(10).name() == "Ca++(aq)"  );
    CHECK( system.species(11).name() == "Mg++(aq)"  );
    CHECK( system.species(12).name() == "CO2(aq)"   );
    CHECK( system.species(13).name() == "HCO3-(aq)" );
    CHECK( system.species(14).name() == "CO3--(aq)" );
    CHECK( system.species(15).name() == "CaCl2(aq)" );
    CHECK( system.species(16).name() == "MgCl2(aq)" );
    CHECK( system.species(17).name() == "SiO2(aq)"  );
    CHECK( system.species(18).name() == "e-(aq)"    );
    CHECK( system.species(19).name() == "CO2(g)"    );
    CHECK( system.species(20).name() == "O2(g)"     );
    CHECK( system.species(21).name() == "H2(g)"     );
    CHECK( system.species(22).name() == "H2O(g)"    );
    CHECK( system.species(23).name() == "CH4(g)"    );
    CHECK( system.species(24).name() == "CO(g)"     );
    CHECK( system.species(25).name() == "NaCl(s)"   );
    CHECK( system.species(26).name() == "CaCO3(s)"  );
    CHECK( system.species(27).name() == "MgCO3(s)"  );
    CHECK( system.species(28).name() == "SiO2(s)"   );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::elements()
    //-------------------------------------------------------------------------
    CHECK( system.elements().size() == 8 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::element(i)
    //-------------------------------------------------------------------------
    Strings symbols = { "H", "O", "Na", "Cl", "C", "Ca", "Mg", "Si" };

    CHECK( contains(symbols, system.element(0).symbol()) );
    CHECK( contains(symbols, system.element(1).symbol()) );
    CHECK( contains(symbols, system.element(2).symbol()) );
    CHECK( contains(symbols, system.element(3).symbol()) );
    CHECK( contains(symbols, system.element(4).symbol()) );
    CHECK( contains(symbols, system.element(5).symbol()) );
    CHECK( contains(symbols, system.element(6).symbol()) );
    CHECK( contains(symbols, system.element(7).symbol()) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::reactions()
    //-------------------------------------------------------------------------
    CHECK( system.reactions().size() == 5 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::reaction(i)
    //-------------------------------------------------------------------------
    CHECK( system.reaction(0).name() == "H2O(aq) = H+(aq) + OH-(aq)" );
    CHECK( system.reaction(1).name() == "H2O(aq) = H2(aq) + 0.5*O2(aq)" );
    CHECK( system.reaction(2).name() == "O2(g) = O2(aq)" );
    CHECK( system.reaction(3).name() == "NaCl(s) = Na+(aq) + Cl-(aq)" );
    CHECK( system.reaction(4).name() == "CaCO3(s)" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::formulaMatrix()
    //-------------------------------------------------------------------------

    const auto Aexpected = detail::assembleFormulaMatrix(system.species(), system.elements());

    CHECK(system.formulaMatrix() == Aexpected);
    CHECK(system.formulaMatrixElements() == Aexpected.topRows(system.elements().size()));
    CHECK(system.formulaMatrixCharge() == Aexpected.row(system.elements().size()));

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::stoichiometricMatrix()
    //-------------------------------------------------------------------------

    const auto Sexpected = detail::assembleStoichiometricMatrix(system.reactions(), system.species());

    CHECK(system.stoichiometricMatrix() == Sexpected);

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::surfaces()
    //-------------------------------------------------------------------------

    auto surfaces = system.surfaces();

    CHECK( surfaces.size() == 3 ); // aqueous-gaseous, aqueous-halite, calcite-implicit phase

    CHECK( surfaces[0].name()         == "AqueousPhase:GaseousPhase" );
    CHECK( surfaces[0].phaseNames()   == Pair<String, String>{ "AqueousPhase", "GaseousPhase" } );
    CHECK( surfaces[0].phaseIndices() == Pair<Index, Index>{ 0, 1 } );

    CHECK( surfaces[1].name()         == "AqueousPhase:Halite" );
    CHECK( surfaces[1].phaseNames()   == Pair<String, String>{ "AqueousPhase", "Halite" } );
    CHECK( surfaces[1].phaseIndices() == Pair<Index, Index>{ 0, 2 } );

    CHECK( surfaces[2].name()         == "Calcite" );
    CHECK( surfaces[2].phaseNames()   == Pair<String, String>{ "Calcite", "Calcite" } );
    CHECK( surfaces[2].phaseIndices() == Pair<Index, Index>{ 3, 3 } );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::surfaces().indexWithPhases()
    //-------------------------------------------------------------------------

    CHECK( system.surfaces().findWithPhases("AqueousPhase", "GaseousPhase") == 0 );
    CHECK( system.surfaces().findWithPhases("AqueousPhase", "Halite") == 1 );
    CHECK( system.surfaces().findWithPhases("Calcite", "Calcite") == 2 );

    CHECK( system.surfaces().findWithPhases("Calcite", "Halite") >= system.surfaces().size() );
    CHECK( system.surfaces().findWithPhases("Calcite", "Quartz") >= system.surfaces().size() );

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: ChemicalSystem::ChemicalSystem(db, phases...)
    //-------------------------------------------------------------------------
    auto db = test::createDatabase();
    system = ChemicalSystem(db,
        AqueousPhase("H2O(aq) H+(aq) OH-(aq) Na+(aq) Cl-(aq)"),
        GaseousPhase("H2O(g)"),
        MineralPhase("NaCl(s)"),
        MineralPhases("CaCO3(s) MgCO3(s) SiO2(s)"));

    CHECK( system.elements().size() == 8 );
    CHECK( system.species().size() == 10 );
    CHECK( system.phases().size() == 6 );
}
