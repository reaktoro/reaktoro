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
#include <Reaktoro/Extensions/Supcrt.hpp>
#include <Reaktoro/Singletons/CriticalProps.hpp>
#include <Reaktoro/Models/ActivityModels.hpp>
#include <Reaktoro/Utils/Material.hpp>
using namespace Reaktoro;

namespace test { auto createChemicalSystem() -> ChemicalSystem; }

auto getSubstanceAmount(Material const& material, ChemicalFormula const& substance)
{
    for(auto [subs, amount] : material.substances())
        if(subs == substance)
            return amount;
    return 0.0;
}

auto getSubstanceMass(Material const& material, ChemicalFormula const& substance)
{
    return getSubstanceAmount(material, substance) * substance.molarMass();
}

auto getSpeciesAmount(Material const& material, String const& species)
{
    const auto idx = material.system().species().index(species);
    for(auto [ispecies, amount] : material.species())
        if(ispecies == idx)
            return amount;
    return 0.0;
}

auto getSpeciesMass(Material const& material, String const& species)
{
    const auto idx = material.system().species().index(species);
    if(idx < material.system().species().size())
        return getSpeciesAmount(material, species) * material.system().species(idx).molarMass();
    else return 0.0;
}

auto getElementAmount(Material const& material, String const& element)
{
    const auto idx = material.system().elements().index(element);
    if(idx < material.system().elements().size())
        return material.elementAmounts()[idx];
    else return 0.0;
}

TEST_CASE("Testing Material class", "[Material]")
{
    const auto system = test::createChemicalSystem();

    Material material(system), othermaterial(system);

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::addSpeciesAmount
    //-------------------------------------------------------------------------
    CHECK( getSpeciesAmount(material, "Ca++(aq)") == Approx(0.0) );

    material.addSpeciesAmount("Ca++(aq)", 1.0, "mol");

    CHECK( getSpeciesAmount(material, "Ca++(aq)") == Approx(1.0) );

    CHECK_THROWS( material.addSpeciesAmount("CH4(aq)", 1.0, "mmol") );  // Species CH4(aq) does not exist in the system

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::addSpeciesMass
    //-------------------------------------------------------------------------
    CHECK( getSpeciesMass(material, "H2O(aq)") == Approx(0.0) );

    material.addSpeciesMass("H2O(aq)", 1.0, "kg");

    CHECK( getSpeciesMass(material, "H2O(aq)") == Approx(1.0) );

    CHECK_THROWS( material.addSpeciesMass("CH3COOH(aq)", 1.0, "mg") );  // Species CH3COOH(aq) does not exist in the system

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::addSubstanceAmount
    //-------------------------------------------------------------------------
    CHECK( getSubstanceAmount(material, "Ca++") == Approx(0.0) );

    material.addSubstanceAmount("Ca++", 1.0, "mol");
    material.addSubstanceAmount("Ca++", 2.0, "mol");

    CHECK( getSubstanceAmount(material, "Ca++") == Approx(3.0) );

    CHECK_THROWS( material.addSubstanceAmount("BrCl2", 1.0, "mol") ); // Substance contains element Br which is not in the system!

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::addSubstanceMass
    //-------------------------------------------------------------------------
    CHECK( getSubstanceMass(material, "SiO2") == Approx(0.0) );

    material.addSubstanceMass("SiO2", 2.0, "kg");

    CHECK( getSubstanceMass(material, "SiO2") == Approx(2.0) );

    CHECK_THROWS( material.addSubstanceMass("ZnCl2", 1.0, "mol") ); // Substance contains element Zn which is not in the system!

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::addMaterialAmount
    //-------------------------------------------------------------------------
    othermaterial = Material(system);
    othermaterial.addMaterialAmount(material, 1.0, "mmol");
    othermaterial.addMaterialAmount(material, 2.0, "mmol");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.003/material.amount() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.003/material.amount() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.003/material.amount() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.003/material.amount() * getSubstanceMass(material, "SiO2") ) );

    ChemicalSystem empty_system;
    Material empty_material(empty_system);

    CHECK_THROWS( othermaterial.addMaterialAmount(empty_material, 1.0, "mmol") ); // Material objects must have the same elements but one is empty!

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::addMaterialMass
    //-------------------------------------------------------------------------
    othermaterial = Material(system);
    othermaterial.addMaterialMass(material, 1.0, "g");
    othermaterial.addMaterialMass(material, 2.0, "g");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.003/material.mass() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.003/material.mass() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.003/material.mass() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.003/material.mass() * getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::add(String, value, unit)
    //-------------------------------------------------------------------------
    othermaterial = Material(system);
    othermaterial.add("Ca++(aq)", 1.0, "mol");
    othermaterial.add("H2O(aq)", 1.0, "kg");
    othermaterial.add("Ca++", 1.0, "mol");
    othermaterial.add("Ca++", 2.0, "mol");
    othermaterial.add("SiO2", 2.0, "kg");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::add(Material, value, unit)
    //-------------------------------------------------------------------------
    othermaterial = Material(system);
    othermaterial.add(material, 1.0, "g");
    othermaterial.add(material, 2.0, "g");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.003/material.mass() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.003/material.mass() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.003/material.mass() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.003/material.mass() * getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::scaleAmount
    //-------------------------------------------------------------------------
    othermaterial = Material(material);
    othermaterial.scaleAmount(7.0, "mmol");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.007/material.amount() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.007/material.amount() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.007/material.amount() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.007/material.amount() * getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::scaleMass
    //-------------------------------------------------------------------------
    othermaterial = Material(material);
    othermaterial.scaleMass(9.0, "g");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.009/material.mass() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.009/material.mass() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.009/material.mass() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.009/material.mass() * getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::scale (using mol based unit)
    //-------------------------------------------------------------------------
    othermaterial = Material(material);
    othermaterial.scaleAmount(11.0, "mmol");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.011/material.amount() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.011/material.amount() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.011/material.amount() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.011/material.amount() * getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::scale (using kg based unit)
    //-------------------------------------------------------------------------
    othermaterial = Material(material);
    othermaterial.scale(6000.0, "mg");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.006/material.mass() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.006/material.mass() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.006/material.mass() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.006/material.mass() * getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::with (using mol based unit)
    //-------------------------------------------------------------------------
    othermaterial = material.with(11.0, "mmol");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.011/material.amount() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.011/material.amount() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.011/material.amount() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.011/material.amount() * getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::with (using kg based unit)
    //-------------------------------------------------------------------------
    othermaterial = material.with(6000.0, "mg");

    CHECK( getSpeciesAmount(othermaterial, "Ca++(aq)") == Approx(0.006/material.mass() * getSpeciesAmount(material, "Ca++(aq)") ) );
    CHECK( getSpeciesMass(othermaterial, "H2O(aq)")    == Approx(0.006/material.mass() * getSpeciesMass(material, "H2O(aq)") ) );
    CHECK( getSubstanceAmount(othermaterial, "Ca++")   == Approx(0.006/material.mass() * getSubstanceAmount(material, "Ca++") ) );
    CHECK( getSubstanceMass(othermaterial, "SiO2")     == Approx(0.006/material.mass() * getSubstanceMass(material, "SiO2") ) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::elementAmounts
    //-------------------------------------------------------------------------
    material = Material(system);
    material.add("H2O"    , 1.0, "mol");
    material.add("CO2(aq)", 2.0, "mol");
    material.add("SiO2"   , 3.0, "mol");

    CHECK( getElementAmount(material, "H" ) == Approx( 2.0) );
    CHECK( getElementAmount(material, "O" ) == Approx(11.0) );
    CHECK( getElementAmount(material, "C" ) == Approx( 2.0) );
    CHECK( getElementAmount(material, "Si") == Approx( 3.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::charge
    //-------------------------------------------------------------------------
    material = Material(system);

    CHECK( material.charge() == Approx(0.0) );

    material.add("Ca++", 1.0, "mol");

    CHECK( material.charge() == Approx(2.0) );

    material.add("Na+(aq)", 1.0, "mol");

    CHECK( material.charge() == Approx(3.0) );

    material.add("Cl-(aq)", 1.0, "mol");

    CHECK( material.charge() == Approx(2.0) );

    material.add("HCO3-", 2.0, "mol");

    CHECK( material.charge() == Approx(0.0) );

    material.add("CO3-2", 2.0, "mol");

    CHECK( material.charge() == Approx(-4.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::componentAmounts
    //-------------------------------------------------------------------------
    const auto E = system.elements().size();

    CHECK( material.componentAmounts().head(E).isApprox(material.elementAmounts()) );
    CHECK( material.componentAmounts()[E] == Approx(material.charge()) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::amount and Material::mass
    //-------------------------------------------------------------------------
    const auto mmCaCO3 = ChemicalFormula("CaCO3").molarMass();
    const auto mmCO2   = ChemicalFormula("CO2").molarMass();
    const auto mmH2O   = ChemicalFormula("H2O").molarMass();
    const auto mmSiO2  = ChemicalFormula("SiO2").molarMass();

    material = Material(system);

    CHECK( material.amount() == Approx(0.0) );
    CHECK( material.mass()   == Approx(0.0) );

    material.add("H2O(aq)", 55.0, "mol");

    CHECK( material.amount() == Approx(55.0) );
    CHECK( material.mass()   == Approx(55.0*mmH2O) );

    material.add("CO2", 10.0, "mol");

    CHECK( material.amount() == Approx(65.0) );
    CHECK( material.mass()   == Approx(55.0*mmH2O + 10.0*mmCO2) );

    material.add("SiO2", 2000.0, "g");

    CHECK( material.amount() == Approx(65.0 + 2.0/mmSiO2) );
    CHECK( material.mass()   == Approx(55.0*mmH2O + 10.0*mmCO2 + 2.0) );

    material.add("CaCO3", 3000.0, "g");

    CHECK( material.amount() == Approx(65.0 + 2.0/mmSiO2 + 3.0/mmCaCO3) );
    CHECK( material.mass()   == Approx(55.0*mmH2O + 10.0*mmCO2 + 2.0 + 3.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::molarMass
    //-------------------------------------------------------------------------
    CHECK( material.molarMass() == Approx(material.mass() / material.amount()) );
}

TEST_CASE("Testing equilibrium capabilities in Material class", "[Material]")
{
    CriticalProps::setMissingAs("He");

    SupcrtDatabase db("supcrtbl");

    AqueousPhase solution(speciate("H O Na Cl Ca C Mg Si"), exclude("organic"));
    solution.setActivityModel(chain(ActivityModelHKF(), ActivityModelDrummond("CO2")));

    GaseousPhase gases("H2O(g) CO2(g) O2(g) CH4(g)");
    gases.setActivityModel(ActivityModelPengRobinson());

    MineralPhases minerals("Halite Calcite Magnesite Dolomite Quartz");

    ChemicalSystem system(db, solution, gases, minerals);

    Material brine(system);
    brine.add("H2O"    , 0.30, "kg");  // using formula
    brine.add("H2O(aq)", 0.70, "kg");  // using species name
    brine.add("NaCl"   , 1.00, "mol"); // using formula
    brine.add("CaCl2"  , 0.10, "mol"); // using formula
    brine.add("MgCl2"  , 0.05, "mol"); // using formula

    Material rock(system);
    rock.add("CaCO3", 1.00, "mol");    // using formula
    rock.add("SiO2" , 1.00, "mol");    // using formula

    Material mix = brine + rock;

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::initialState
    //-------------------------------------------------------------------------
    ChemicalState state0(system);

    state0 = mix.initialState(300.0, 10.0e5); // at 300 K and 10 bar

    CHECK( state0.temperature() == Approx(300.0) );
    CHECK( state0.pressure() == Approx(10e5) );

    CHECK( state0.speciesMass("H2O(aq)")     == Approx(1.00) ); // because of brine.add("H2O", 0.30, "kg") and brine.add("H2O(aq)", 0.70, "kg")
    CHECK( state0.speciesAmount("NaCl(aq)")  == Approx(1.00) ); // because of brine.add("NaCl" , 1.00, "mol")
    CHECK( state0.speciesAmount("CaCl2(aq)") == Approx(0.10) ); // because of brine.add("CaCl2" , 0.10, "mol")
    CHECK( state0.speciesAmount("Calcite")   == Approx(1.00) ); // because of rock.add("CaCO3", 1.00, "mol")
    CHECK( state0.speciesAmount("Quartz")    == Approx(1.00) ); // because of rock.add("SiO2", 1.00, "mol")

    //-------------------------------------------------------------------------
    // TESTING METHOD: Material::equilibrate
    //-------------------------------------------------------------------------
    mix = brine.with(1.0, "kg") + rock.with(1.0, "kg");

    const ArrayXr bmix = mix.componentAmounts(); // the amounts of elements/charge in mix material

    ChemicalState state(system);

    state = mix.equilibrate();

    CHECK( mix.result().succeeded() );
    CHECK( mix.result().iterations() == 76 );

    CHECK( state.temperature() == Approx(25.0 + 273.15) );
    CHECK( state.pressure() == Approx(1.0 * 1e5) );
    CHECK( state.componentAmounts().isApprox(bmix) ); // ensure computed state has element/charge amounts equal to those in mix material

    state = mix.equilibrate(60.0, "celsius", 10.0, "bar");

    CHECK( mix.result().succeeded() );
    CHECK( mix.result().iterations() == 90 );

    CHECK( state.temperature() == Approx(60.0 + 273.15) );
    CHECK( state.pressure() == Approx(10.0 * 1e5) );
    CHECK( state.componentAmounts().isApprox(bmix) ); // ensure computed state has element/charge amounts equal to those in mix material
}
