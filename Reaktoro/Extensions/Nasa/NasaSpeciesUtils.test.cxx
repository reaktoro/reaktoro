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
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Extensions/Nasa/NasaSpecies.hpp>
#include <Reaktoro/Extensions/Nasa/NasaSpeciesUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing NasaSpeciesUtils module", "[NasaSpeciesUtils]")
{
    using namespace NasaUtils;

    //======================================================================
    // Testing method NasaUtils::correctElementSymbol
    //======================================================================

    CHECK( correctElementSymbol("AG") == "Ag" );
    CHECK( correctElementSymbol("AL") == "Al" );
    CHECK( correctElementSymbol("AR") == "Ar" );
    CHECK( correctElementSymbol("B")  == "B"  );
    CHECK( correctElementSymbol("BA") == "Ba" );
    CHECK( correctElementSymbol("BE") == "Be" );
    CHECK( correctElementSymbol("BR") == "Br" );
    CHECK( correctElementSymbol("C")  == "C"  );
    CHECK( correctElementSymbol("CA") == "Ca" );
    CHECK( correctElementSymbol("CD") == "Cd" );
    CHECK( correctElementSymbol("CL") == "Cl" );
    CHECK( correctElementSymbol("CO") == "Co" );
    CHECK( correctElementSymbol("CR") == "Cr" );
    CHECK( correctElementSymbol("CS") == "Cs" );
    CHECK( correctElementSymbol("CU") == "Cu" );
    CHECK( correctElementSymbol("D")  == "D"  );
    CHECK( correctElementSymbol("E")  == "E"  );
    CHECK( correctElementSymbol("F")  == "F"  );
    CHECK( correctElementSymbol("FE") == "Fe" );
    CHECK( correctElementSymbol("GA") == "Ga" );
    CHECK( correctElementSymbol("GE") == "Ge" );
    CHECK( correctElementSymbol("H")  == "H"  );
    CHECK( correctElementSymbol("HE") == "He" );
    CHECK( correctElementSymbol("HG") == "Hg" );
    CHECK( correctElementSymbol("I")  == "I"  );
    CHECK( correctElementSymbol("IN") == "In" );
    CHECK( correctElementSymbol("K")  == "K"  );
    CHECK( correctElementSymbol("KR") == "Kr" );
    CHECK( correctElementSymbol("LI") == "Li" );
    CHECK( correctElementSymbol("MG") == "Mg" );
    CHECK( correctElementSymbol("MN") == "Mn" );
    CHECK( correctElementSymbol("MO") == "Mo" );
    CHECK( correctElementSymbol("N")  == "N"  );
    CHECK( correctElementSymbol("NA") == "Na" );
    CHECK( correctElementSymbol("NB") == "Nb" );
    CHECK( correctElementSymbol("NE") == "Ne" );
    CHECK( correctElementSymbol("NI") == "Ni" );
    CHECK( correctElementSymbol("O")  == "O"  );
    CHECK( correctElementSymbol("P")  == "P"  );
    CHECK( correctElementSymbol("PB") == "Pb" );
    CHECK( correctElementSymbol("RB") == "Rb" );
    CHECK( correctElementSymbol("RN") == "Rn" );
    CHECK( correctElementSymbol("S")  == "S"  );
    CHECK( correctElementSymbol("SC") == "Sc" );
    CHECK( correctElementSymbol("SI") == "Si" );
    CHECK( correctElementSymbol("SN") == "Sn" );
    CHECK( correctElementSymbol("SR") == "Sr" );
    CHECK( correctElementSymbol("TA") == "Ta" );
    CHECK( correctElementSymbol("TH") == "Th" );
    CHECK( correctElementSymbol("TI") == "Ti" );
    CHECK( correctElementSymbol("U")  == "U"  );
    CHECK( correctElementSymbol("V")  == "V"  );
    CHECK( correctElementSymbol("W")  == "W"  );
    CHECK( correctElementSymbol("XE") == "Xe" );
    CHECK( correctElementSymbol("ZN") == "Zn" );
    CHECK( correctElementSymbol("ZR") == "Zr" );

    //======================================================================
    // Testing method NasaUtils::createElement
    //======================================================================

    CHECK_NOTHROW( createElement("AG") );
    CHECK_NOTHROW( createElement("AL") );
    CHECK_NOTHROW( createElement("AR") );
    CHECK_NOTHROW( createElement("B")  );
    CHECK_NOTHROW( createElement("BA") );
    CHECK_NOTHROW( createElement("BE") );
    CHECK_NOTHROW( createElement("BR") );
    CHECK_NOTHROW( createElement("C")  );
    CHECK_NOTHROW( createElement("CA") );
    CHECK_NOTHROW( createElement("CD") );
    CHECK_NOTHROW( createElement("CL") );
    CHECK_NOTHROW( createElement("CO") );
    CHECK_NOTHROW( createElement("CR") );
    CHECK_NOTHROW( createElement("CS") );
    CHECK_NOTHROW( createElement("CU") );
    CHECK_NOTHROW( createElement("D")  );
    CHECK_NOTHROW( createElement("F")  );
    CHECK_NOTHROW( createElement("FE") );
    CHECK_NOTHROW( createElement("GA") );
    CHECK_NOTHROW( createElement("GE") );
    CHECK_NOTHROW( createElement("H")  );
    CHECK_NOTHROW( createElement("HE") );
    CHECK_NOTHROW( createElement("HG") );
    CHECK_NOTHROW( createElement("I")  );
    CHECK_NOTHROW( createElement("IN") );
    CHECK_NOTHROW( createElement("K")  );
    CHECK_NOTHROW( createElement("KR") );
    CHECK_NOTHROW( createElement("LI") );
    CHECK_NOTHROW( createElement("MG") );
    CHECK_NOTHROW( createElement("MN") );
    CHECK_NOTHROW( createElement("MO") );
    CHECK_NOTHROW( createElement("N")  );
    CHECK_NOTHROW( createElement("NA") );
    CHECK_NOTHROW( createElement("NB") );
    CHECK_NOTHROW( createElement("NE") );
    CHECK_NOTHROW( createElement("NI") );
    CHECK_NOTHROW( createElement("O")  );
    CHECK_NOTHROW( createElement("P")  );
    CHECK_NOTHROW( createElement("PB") );
    CHECK_NOTHROW( createElement("RB") );
    CHECK_NOTHROW( createElement("RN") );
    CHECK_NOTHROW( createElement("S")  );
    CHECK_NOTHROW( createElement("SC") );
    CHECK_NOTHROW( createElement("SI") );
    CHECK_NOTHROW( createElement("SN") );
    CHECK_NOTHROW( createElement("SR") );
    CHECK_NOTHROW( createElement("TA") );
    CHECK_NOTHROW( createElement("TH") );
    CHECK_NOTHROW( createElement("TI") );
    CHECK_NOTHROW( createElement("U")  );
    CHECK_NOTHROW( createElement("V")  );
    CHECK_NOTHROW( createElement("W")  );
    CHECK_NOTHROW( createElement("XE") );
    CHECK_NOTHROW( createElement("ZN") );
    CHECK_NOTHROW( createElement("ZR") );

    CHECK_THROWS( createElement("E")  ); // charge element symbol in NASA database - not needed in Reaktoro as we use Species::withCharge method

    //======================================================================
    // Testing method NasaUtils::createElements
    //======================================================================

    NasaSpecies species01;
    species01.name = "MgO";
    species01.formula = {{"MG", 1.0}, {"O", 1.0}};
    species01.molarmass = 40.3044000;
    species01.aggregatestate = NasaAggregateState::Gas;
    species01.type = NasaSpeciesType::Product;

    NasaSpecies species02;
    species02.name = "Mg(cr)";
    species02.formula = {{"MG", 1.0}};
    species02.molarmass = 24.3050000;
    species02.aggregatestate = NasaAggregateState::Condensed;
    species02.type = NasaSpeciesType::Product;

    NasaSpecies species03;
    species03.name = "Jet-A(g)";
    species03.formula = {{"C", 12.0}, {"H", 23.0}};
    species03.molarmass = 167.3110200;
    species03.aggregatestate = NasaAggregateState::Gas;
    species03.type = NasaSpeciesType::Product;

    NasaSpecies species04;
    species04.name = "NaCN(II)";
    species04.formula = {{"NA", 1.0}, {"C", 1.0}, {"N", 1.0}};
    species04.molarmass = 49.0071700;
    species04.aggregatestate = NasaAggregateState::Condensed;
    species04.type = NasaSpeciesType::Product;

    NasaSpecies species05;
    species05.name = "C2H2(L),acetyle";
    species05.formula = {{"C", 2.0}, {"H", 2.0}};
    species05.molarmass = 26.0372800;
    species05.aggregatestate = NasaAggregateState::Condensed;
    species05.type = NasaSpeciesType::Product;

    NasaSpecies species06;
    species06.name = "N2O4(L)";
    species06.formula = {{"N", 2.0}, {"O", 4.0}};
    species06.molarmass = 92.0110000;
    species06.aggregatestate = NasaAggregateState::Condensed;
    species06.type = NasaSpeciesType::Product;

    NasaSpecies species07;
    species07.name = "C6D5,phenyl";
    species07.formula = {{"C", 6.0}, {"D", 5.0}};
    species07.molarmass = 82.1347100;
    species07.aggregatestate = NasaAggregateState::Gas;
    species07.type = NasaSpeciesType::Reactant;

    NasaSpecies species08;
    species08.name = "CHFCLBr";
    species08.formula = {{"C", 1.0}, {"H", 1.0}, {"F", 1.0}, {"CL", 1.0}, {"BR", 1.0}};
    species08.molarmass = 147.3740432;
    species08.aggregatestate = NasaAggregateState::Gas;
    species08.type = NasaSpeciesType::Reactant;

    NasaSpecies species09;
    species09.name = "SrOH+";
    species09.formula = {{"SR", 1.0}, {"O", 1.0}, {"H", 1.0}, {"E", -1.0}};
    species09.molarmass = 104.6267914;
    species09.aggregatestate = NasaAggregateState::Gas;
    species09.type = NasaSpeciesType::Reactant;

    NasaSpecies species10;
    species10.name = "Be++";
    species10.formula = {{"BE", 1.0}, {"E", -2.0}};
    species10.molarmass = 9.0110848;
    species10.aggregatestate = NasaAggregateState::Gas;
    species10.type = NasaSpeciesType::Reactant;

    NasaSpecies species11;
    species11.name = "Br-";
    species11.formula = {{"BR", 1.0}, {"E", 1.0}};
    species11.molarmass = 79.9045486;
    species11.aggregatestate = NasaAggregateState::Gas;
    species11.type = NasaSpeciesType::Reactant;

    const ElementalComposition elements01 = createElements(species01);
    const ElementalComposition elements02 = createElements(species02);
    const ElementalComposition elements03 = createElements(species03);
    const ElementalComposition elements04 = createElements(species04);
    const ElementalComposition elements05 = createElements(species05);
    const ElementalComposition elements06 = createElements(species06);
    const ElementalComposition elements07 = createElements(species07);
    const ElementalComposition elements08 = createElements(species08);
    const ElementalComposition elements09 = createElements(species09);
    const ElementalComposition elements10 = createElements(species10);
    const ElementalComposition elements11 = createElements(species11);

    CHECK( elements01.size() == 2 );
    CHECK( elements01.coefficient("Mg") == 1.0 );
    CHECK( elements01.coefficient("O")  == 1.0 );

    CHECK( elements02.size() == 1 );
    CHECK( elements02.coefficient("Mg") == 1.0 );

    CHECK( elements03.size() == 2 );
    CHECK( elements03.coefficient("C") == 12.0 );
    CHECK( elements03.coefficient("H") == 23.0 );

    CHECK( elements04.size() == 3 );
    CHECK( elements04.coefficient("Na") == 1.0 );
    CHECK( elements04.coefficient("C")  == 1.0 );
    CHECK( elements04.coefficient("N")  == 1.0 );

    CHECK( elements05.size() == 2 );
    CHECK( elements05.coefficient("C") == 2.0 );
    CHECK( elements05.coefficient("H") == 2.0 );

    CHECK( elements06.size() == 2 );
    CHECK( elements06.coefficient("N") == 2.0 );
    CHECK( elements06.coefficient("O") == 4.0 );

    CHECK( elements07.size() == 2 );
    CHECK( elements07.coefficient("C") == 6.0 );
    CHECK( elements07.coefficient("D") == 5.0 );

    CHECK( elements08.size() == 5 );
    CHECK( elements08.coefficient("C")  == 1.0 );
    CHECK( elements08.coefficient("H")  == 1.0 );
    CHECK( elements08.coefficient("F")  == 1.0 );
    CHECK( elements08.coefficient("Cl") == 1.0 );
    CHECK( elements08.coefficient("Br") == 1.0 );

    CHECK( elements09.size() == 3 );
    CHECK( elements09.coefficient("Sr") == 1.0  );
    CHECK( elements09.coefficient("H")  == 1.0  );
    CHECK( elements09.coefficient("O")  == 1.0  );

    CHECK( elements10.size() == 1 );
    CHECK( elements10.coefficient("Be") == 1.0  );

    //======================================================================
    // Testing method NasaUtils::charge
    //======================================================================

    CHECK( NasaUtils::charge(species01) ==  0.0 );
    CHECK( NasaUtils::charge(species02) ==  0.0 );
    CHECK( NasaUtils::charge(species03) ==  0.0 );
    CHECK( NasaUtils::charge(species04) ==  0.0 );
    CHECK( NasaUtils::charge(species05) ==  0.0 );
    CHECK( NasaUtils::charge(species06) ==  0.0 );
    CHECK( NasaUtils::charge(species07) ==  0.0 );
    CHECK( NasaUtils::charge(species08) ==  0.0 );
    CHECK( NasaUtils::charge(species09) ==  1.0 );
    CHECK( NasaUtils::charge(species10) ==  2.0 );
    CHECK( NasaUtils::charge(species11) == -1.0 );

    //======================================================================
    // Testing method NasaUtils::aggregateState
    //======================================================================

    CHECK( NasaUtils::aggregateState(species01) == AggregateState::Gas );
    CHECK( NasaUtils::aggregateState(species02) == AggregateState::CondensedPhase );
    CHECK( NasaUtils::aggregateState(species03) == AggregateState::Gas );
    CHECK( NasaUtils::aggregateState(species04) == AggregateState::CondensedPhase );
    CHECK( NasaUtils::aggregateState(species05) == AggregateState::CondensedPhase );
    CHECK( NasaUtils::aggregateState(species06) == AggregateState::CondensedPhase );
    CHECK( NasaUtils::aggregateState(species07) == AggregateState::Gas );
    CHECK( NasaUtils::aggregateState(species08) == AggregateState::Gas );
    CHECK( NasaUtils::aggregateState(species09) == AggregateState::Gas );
    CHECK( NasaUtils::aggregateState(species10) == AggregateState::Gas );
    CHECK( NasaUtils::aggregateState(species11) == AggregateState::Gas );

    //======================================================================
    // Testing method NasaUtils::tags
    //======================================================================

    CHECK( NasaUtils::tags(species01) == Strings{ "product" } );
    CHECK( NasaUtils::tags(species02) == Strings{ "product" } );
    CHECK( NasaUtils::tags(species03) == Strings{ "product" } );
    CHECK( NasaUtils::tags(species04) == Strings{ "product" } );
    CHECK( NasaUtils::tags(species05) == Strings{ "product" } );
    CHECK( NasaUtils::tags(species06) == Strings{ "product" } );
    CHECK( NasaUtils::tags(species07) == Strings{ "reactant" } );
    CHECK( NasaUtils::tags(species08) == Strings{ "reactant" } );
    CHECK( NasaUtils::tags(species09) == Strings{ "reactant" } );
    CHECK( NasaUtils::tags(species10) == Strings{ "reactant" } );
    CHECK( NasaUtils::tags(species11) == Strings{ "reactant" } );

    //======================================================================
    // Testing method NasaUtils::convertSpecies
    //======================================================================

    Species s01 = NasaUtils::convertSpecies(species01);
    Species s02 = NasaUtils::convertSpecies(species02);
    Species s03 = NasaUtils::convertSpecies(species03);
    Species s04 = NasaUtils::convertSpecies(species04);
    Species s05 = NasaUtils::convertSpecies(species05);
    Species s06 = NasaUtils::convertSpecies(species06);
    Species s07 = NasaUtils::convertSpecies(species07);
    Species s08 = NasaUtils::convertSpecies(species08);
    Species s09 = NasaUtils::convertSpecies(species09);
    Species s10 = NasaUtils::convertSpecies(species10);
    Species s11 = NasaUtils::convertSpecies(species11);

    CHECK( s01.name() == species01.name );
    CHECK( s02.name() == species02.name );
    CHECK( s03.name() == species03.name );
    CHECK( s04.name() == species04.name );
    CHECK( s05.name() == species05.name );
    CHECK( s06.name() == species06.name );
    CHECK( s07.name() == species07.name );
    CHECK( s08.name() == species08.name );
    CHECK( s09.name() == species09.name );
    CHECK( s10.name() == species10.name );
    CHECK( s11.name() == species11.name );

    CHECK( s01.formula().str() == species01.name );
    CHECK( s02.formula().str() == species02.name );
    CHECK( s03.formula().str() == species03.name );
    CHECK( s04.formula().str() == species04.name );
    CHECK( s05.formula().str() == species05.name );
    CHECK( s06.formula().str() == species06.name );
    CHECK( s07.formula().str() == species07.name );
    CHECK( s08.formula().str() == species08.name );
    CHECK( s09.formula().str() == species09.name );
    CHECK( s10.formula().str() == species10.name );
    CHECK( s11.formula().str() == species11.name );

    CHECK( s01.elements() == elements01 );
    CHECK( s02.elements() == elements02 );
    CHECK( s03.elements() == elements03 );
    CHECK( s04.elements() == elements04 );
    CHECK( s05.elements() == elements05 );
    CHECK( s06.elements() == elements06 );
    CHECK( s07.elements() == elements07 );
    CHECK( s08.elements() == elements08 );
    CHECK( s09.elements() == elements09 );
    CHECK( s10.elements() == elements10 );
    CHECK( s11.elements() == elements11 );

    CHECK( s01.charge() == NasaUtils::charge(species01) );
    CHECK( s02.charge() == NasaUtils::charge(species02) );
    CHECK( s03.charge() == NasaUtils::charge(species03) );
    CHECK( s04.charge() == NasaUtils::charge(species04) );
    CHECK( s05.charge() == NasaUtils::charge(species05) );
    CHECK( s06.charge() == NasaUtils::charge(species06) );
    CHECK( s07.charge() == NasaUtils::charge(species07) );
    CHECK( s08.charge() == NasaUtils::charge(species08) );
    CHECK( s09.charge() == NasaUtils::charge(species09) );
    CHECK( s10.charge() == NasaUtils::charge(species10) );
    CHECK( s11.charge() == NasaUtils::charge(species11) );

    CHECK( s01.aggregateState() == NasaUtils::aggregateState(species01) );
    CHECK( s02.aggregateState() == NasaUtils::aggregateState(species02) );
    CHECK( s03.aggregateState() == NasaUtils::aggregateState(species03) );
    CHECK( s04.aggregateState() == NasaUtils::aggregateState(species04) );
    CHECK( s05.aggregateState() == NasaUtils::aggregateState(species05) );
    CHECK( s06.aggregateState() == NasaUtils::aggregateState(species06) );
    CHECK( s07.aggregateState() == NasaUtils::aggregateState(species07) );
    CHECK( s08.aggregateState() == NasaUtils::aggregateState(species08) );
    CHECK( s09.aggregateState() == NasaUtils::aggregateState(species09) );
    CHECK( s10.aggregateState() == NasaUtils::aggregateState(species10) );
    CHECK( s11.aggregateState() == NasaUtils::aggregateState(species11) );

    CHECK( s01.molarMass() == Approx(species01.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s02.molarMass() == Approx(species02.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s03.molarMass() == Approx(species03.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s04.molarMass() == Approx(species04.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s05.molarMass() == Approx(species05.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s06.molarMass() == Approx(species06.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s07.molarMass() == Approx(species07.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s08.molarMass() == Approx(species08.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s09.molarMass() == Approx(species09.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s10.molarMass() == Approx(species10.molarmass * 1e-3).epsilon(1e-3) );
    CHECK( s11.molarMass() == Approx(species11.molarmass * 1e-3).epsilon(1e-3) );

    CHECK( s01.standardThermoModel() );
    CHECK( s02.standardThermoModel() );
    CHECK( s03.standardThermoModel() );
    CHECK( s04.standardThermoModel() );
    CHECK( s05.standardThermoModel() );
    CHECK( s06.standardThermoModel() );
    CHECK( s07.standardThermoModel() );
    CHECK( s08.standardThermoModel() );
    CHECK( s09.standardThermoModel() );
    CHECK( s10.standardThermoModel() );
    CHECK( s11.standardThermoModel() );
}

