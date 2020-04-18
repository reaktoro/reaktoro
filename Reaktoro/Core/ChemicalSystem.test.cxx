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
#include <Reaktoro/Core/ChemicalSystem.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ChemicalSystem class", "[ChemicalSystem]")
{
    // Create the Phase objects in the ChemicalSystem
    const Vec<Phase> phases =
    {
        Phase()
            .withName("AqueousSolution")
            .withSpecies(SpeciesList("H2O(aq) H+(aq) OH-(aq) H2(aq) O2(aq) Na+(aq) Cl-(aq) NaCl(aq) HCO3-(aq) CO2(aq) CO3--(aq)"))
            .withStateOfMatter(StateOfMatter::Liquid),
        Phase()
            .withName("GaseousSolution")
            .withSpecies(SpeciesList("H2O(g) CO2(g) H2(g) O2(g)"))
            .withStateOfMatter(StateOfMatter::Gas),
        Phase()
            .withName("Halite")
            .withSpecies(SpeciesList("NaCl(s)"))
            .withStateOfMatter(StateOfMatter::Solid),
        Phase()
            .withName("Calcite")
            .withSpecies(SpeciesList("CaCO3(s)"))
            .withStateOfMatter(StateOfMatter::Solid),
        Phase()
            .withName("Quartz")
            .withSpecies(SpeciesList("SiO2(s)"))
            .withStateOfMatter(StateOfMatter::Solid)
    };

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTOR: ChemicalSystem::ChemicalSystem(phases)
    //-------------------------------------------------------------------------
    ChemicalSystem system(phases);

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::phases()
    //-------------------------------------------------------------------------
    REQUIRE( system.phases().size() == 5 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::phase(i)
    //-------------------------------------------------------------------------
    REQUIRE( system.phase(0).name() == "AqueousSolution" );
    REQUIRE( system.phase(1).name() == "GaseousSolution" );
    REQUIRE( system.phase(2).name() == "Halite"          );
    REQUIRE( system.phase(3).name() == "Calcite"         );
    REQUIRE( system.phase(4).name() == "Quartz"          );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::species()
    //-------------------------------------------------------------------------
    REQUIRE( system.species().size() == 18 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::species(i)
    //-------------------------------------------------------------------------
    REQUIRE( system.species(0).name()  == "H2O(aq)"   );
    REQUIRE( system.species(1).name()  == "H+(aq)"    );
    REQUIRE( system.species(2).name()  == "OH-(aq)"   );
    REQUIRE( system.species(3).name()  == "H2(aq)"    );
    REQUIRE( system.species(4).name()  == "O2(aq)"    );
    REQUIRE( system.species(5).name()  == "Na+(aq)"   );
    REQUIRE( system.species(6).name()  == "Cl-(aq)"   );
    REQUIRE( system.species(7).name()  == "NaCl(aq)"  );
    REQUIRE( system.species(8).name()  == "HCO3-(aq)" );
    REQUIRE( system.species(9).name()  == "CO2(aq)"   );
    REQUIRE( system.species(10).name() == "CO3--(aq)" );
    REQUIRE( system.species(11).name() == "H2O(g)"    );
    REQUIRE( system.species(12).name() == "CO2(g)"    );
    REQUIRE( system.species(13).name() == "H2(g)"     );
    REQUIRE( system.species(14).name() == "O2(g)"     );
    REQUIRE( system.species(15).name() == "NaCl(s)"   );
    REQUIRE( system.species(16).name() == "CaCO3(s)"  );
    REQUIRE( system.species(17).name() == "SiO2(s)"   );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::elements()
    //-------------------------------------------------------------------------
    REQUIRE( system.elements().size() == 7 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::element(i)
    //-------------------------------------------------------------------------
    Strings symbols = { "H", "O", "Na", "Cl", "C", "Ca", "Si" };

    REQUIRE( contains(symbols, system.element(0).symbol()) );
    REQUIRE( contains(symbols, system.element(1).symbol()) );
    REQUIRE( contains(symbols, system.element(2).symbol()) );
    REQUIRE( contains(symbols, system.element(3).symbol()) );
    REQUIRE( contains(symbols, system.element(4).symbol()) );
    REQUIRE( contains(symbols, system.element(5).symbol()) );
    REQUIRE( contains(symbols, system.element(6).symbol()) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: ChemicalSystem::formulaMatrix()
    //-------------------------------------------------------------------------
    // The transpose of the formula matrix of the chemical system.
    const MatrixXd Atr{
       // H   C   O  Na  Si  Cl  Ca   Z     (elements sorted in ascending order of atomic weights)
        { 2,  0,  1,  0,  0,  0,  0,  0 },  // H2O(aq)
        { 1,  0,  0,  0,  0,  0,  0,  1 },  // H+(aq)
        { 1,  0,  1,  0,  0,  0,  0, -1 },  // OH-(aq)
        { 2,  0,  0,  0,  0,  0,  0,  0 },  // H2(aq)
        { 0,  0,  2,  0,  0,  0,  0,  0 },  // O2(aq)
        { 0,  0,  0,  1,  0,  0,  0,  1 },  // Na+(aq)
        { 0,  0,  0,  0,  0,  1,  0, -1 },  // Cl-(aq)
        { 0,  0,  0,  1,  0,  1,  0,  0 },  // NaCl(aq)
        { 1,  1,  3,  0,  0,  0,  0, -1 },  // HCO3-(aq)
        { 0,  1,  2,  0,  0,  0,  0,  0 },  // CO2(aq)
        { 0,  1,  3,  0,  0,  0,  0, -2 },  // CO3--(aq)
        { 2,  0,  1,  0,  0,  0,  0,  0 },  // H2O(g)
        { 0,  1,  2,  0,  0,  0,  0,  0 },  // CO2(g)
        { 2,  0,  0,  0,  0,  0,  0,  0 },  // H2(g)
        { 0,  0,  2,  0,  0,  0,  0,  0 },  // O2(g)
        { 0,  0,  0,  1,  0,  1,  0,  0 },  // NaCl(s)
        { 0,  1,  3,  0,  0,  0,  1,  0 },  // CaCO3(s)
        { 0,  0,  2,  0,  1,  0,  0,  0 }}; // SiO2(s)

    REQUIRE( Atr.transpose() == system.formulaMatrix() );
}