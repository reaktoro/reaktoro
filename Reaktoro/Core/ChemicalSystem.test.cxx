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
using namespace Reaktoro;

namespace test {

/// Return a mock Database object for test reasons (imported).
extern auto createDatabase() -> Database;

/// Return mock activity properties for an aqueous solution.
auto activityPropsAqueous(ActivityPropsRef props, ActivityArgs args)
{
    const auto [T, P, x, extra] = args;
    props.Vex  = 0.1 * log(P/T);
    props.VexT = 0.2 * log(P/T);
    props.VexP = 0.3 * log(P/T);
    props.Gex  = 0.4 * log(P/T);
    props.Hex  = 0.5 * log(P/T);
    props.Cpex = 0.6 * log(P/T);
    props.Cvex = 0.7 * log(P/T);
    props.ln_g = 0.8 * x;
    props.ln_a = 0.9 * x;
};

/// Return mock activity properties for a gaseous solution.
auto activityPropsGaseous(ActivityPropsRef props, ActivityArgs args)
{
    const auto [T, P, x, extra] = args;
    props.Vex  = 1.0 * log(P/T);
    props.VexT = 2.0 * log(P/T);
    props.VexP = 3.0 * log(P/T);
    props.Gex  = 4.0 * log(P/T);
    props.Hex  = 5.0 * log(P/T);
    props.Cpex = 6.0 * log(P/T);
    props.Cvex = 7.0 * log(P/T);
    props.ln_g = 8.0 * x;
    props.ln_a = 9.0 * x;
};

/// Return mock activity properties for a solid solution.
auto activityPropsSolid(ActivityPropsRef props, ActivityArgs args)
{
    const auto [T, P, x, extra] = args;
    props.Vex  = 1.1 * log(P/T);
    props.VexT = 2.1 * log(P/T);
    props.VexP = 3.1 * log(P/T);
    props.Gex  = 4.1 * log(P/T);
    props.Hex  = 5.1 * log(P/T);
    props.Cpex = 6.1 * log(P/T);
    props.Cvex = 7.1 * log(P/T);
    props.ln_g = 8.1 * x;
    props.ln_a = 9.1 * x;
};

/// Return a mock ChemicalSystem object for test reasons.
auto createChemicalSystem() -> ChemicalSystem
{
    // Create the Database object for the ChemicalSystem
    Database db = test::createDatabase();

    // Create the Phase objects for the ChemicalSystem
    const Vec<Phase> phases =
    {
        Phase()
            .withName("AqueousPhase")
            .withSpecies(db.speciesWithAggregateState(AggregateState::Aqueous))
            .withStateOfMatter(StateOfMatter::Liquid)
            .withActivityPropsFn(activityPropsAqueous)
            .withIdealActivityPropsFn(activityPropsAqueous),
        Phase()
            .withName("GaseousPhase")
            .withSpecies(db.speciesWithAggregateState(AggregateState::Gas))
            .withStateOfMatter(StateOfMatter::Gas)
            .withActivityPropsFn(activityPropsGaseous)
            .withIdealActivityPropsFn(activityPropsGaseous),
        Phase()
            .withName("Halite")
            .withSpecies({ db.species().get("NaCl(s)") })
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activityPropsSolid)
            .withIdealActivityPropsFn(activityPropsSolid),
        Phase()
            .withName("Calcite")
            .withSpecies({ db.species().get("CaCO3(s)") })
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activityPropsSolid)
            .withIdealActivityPropsFn(activityPropsSolid),
        Phase()
            .withName("Magnesite")
            .withSpecies({ db.species().get("MgCO3(s)") })
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activityPropsSolid)
            .withIdealActivityPropsFn(activityPropsSolid),
        Phase()
            .withName("Quartz")
            .withSpecies({ db.species().get("SiO2(s)") })
            .withStateOfMatter(StateOfMatter::Solid)
            .withActivityPropsFn(activityPropsSolid)
            .withIdealActivityPropsFn(activityPropsSolid),
    };

    return ChemicalSystem(db, phases);
}

} // namespace test

TEST_CASE("Testing ChemicalSystem class", "[ChemicalSystem]")
{
    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: ChemicalSystem::ChemicalSystem(phases)
    //-------------------------------------------------------------------------
    ChemicalSystem system = test::createChemicalSystem();

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
    // TESTING METHOD: ChemicalSystem::formulaMatrix()
    //-------------------------------------------------------------------------
    // The transpose of the formula matrix of the chemical system.
    const MatrixXd Atr{
       // H   C   O  Na  Mg, Si  Cl  Ca   Z     (elements sorted in ascending order of atomic weights)
        { 2,  0,  1,  0,  0,  0,  0,  0,  0 }, // H2O
        { 1,  0,  0,  0,  0,  0,  0,  0,  1 }, // H+
        { 1,  0,  1,  0,  0,  0,  0,  0, -1 }, // OH-
        { 2,  0,  0,  0,  0,  0,  0,  0,  0 }, // H2
        { 0,  0,  2,  0,  0,  0,  0,  0,  0 }, // O2
        { 0,  0,  0,  1,  0,  0,  0,  0,  1 }, // Na+
        { 0,  0,  0,  0,  0,  0,  1,  0, -1 }, // Cl-
        { 0,  0,  0,  1,  0,  0,  1,  0,  0 }, // NaCl
        { 1,  0,  0,  0,  0,  0,  1,  0,  0 }, // HCl
        { 1,  0,  1,  1,  0,  0,  0,  0,  0 }, // NaOH
        { 0,  0,  0,  0,  0,  0,  0,  1,  2 }, // Ca++
        { 0,  0,  0,  0,  1,  0,  0,  0,  2 }, // Mg++
        { 0,  1,  2,  0,  0,  0,  0,  0,  0 }, // CO2
        { 1,  1,  3,  0,  0,  0,  0,  0, -1 }, // HCO3-
        { 0,  1,  3,  0,  0,  0,  0,  0, -2 }, // CO3--
        { 0,  0,  0,  0,  0,  0,  2,  1,  0 }, // CaCl2
        { 0,  0,  0,  0,  1,  0,  2,  0,  0 }, // MgCl2
        { 0,  0,  2,  0,  0,  1,  0,  0,  0 }, // SiO2
        { 0,  0,  0,  0,  0,  0,  0,  0, -1 }, // e-
        { 0,  1,  2,  0,  0,  0,  0,  0,  0 }, // CO2(g)
        { 0,  0,  2,  0,  0,  0,  0,  0,  0 }, // O2(g)
        { 2,  0,  0,  0,  0,  0,  0,  0,  0 }, // H2(g)
        { 2,  0,  1,  0,  0,  0,  0,  0,  0 }, // H2O(g)
        { 4,  1,  0,  0,  0,  0,  0,  0,  0 }, // CH4(g)
        { 0,  1,  1,  0,  0,  0,  0,  0,  0 }, // CO(g)
        { 0,  0,  0,  1,  0,  0,  1,  0,  0 }, // NaCl(s)
        { 0,  1,  3,  0,  0,  0,  0,  1,  0 }, // CaCO3(s)
        { 0,  1,  3,  0,  1,  0,  0,  0,  0 }, // MgCO3(s)
        { 0,  0,  2,  0,  0,  1,  0,  0,  0 }  // SiO2(s)
    };

    CHECK( Atr.transpose() == system.formulaMatrix() );
}
