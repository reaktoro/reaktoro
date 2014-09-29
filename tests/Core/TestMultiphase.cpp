// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "TestPhase.hpp"

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

namespace Reaktor {
namespace {

// The indices of the species
const unsigned iH2O  = 0;
const unsigned iHp   = 1;
const unsigned iOHm  = 2;
const unsigned iCO2g = 3;
const unsigned iH2Og = 4;

// The indices of the elements
const unsigned iC    = 0;
const unsigned iH    = 1;
const unsigned iO    = 2;

auto createMultiphase() -> Multiphase
{
    ThermoScalar thermo_property(1.0, 2.0, 3.0);
    ThermoScalarFunction thermo_property_fn = [=](double,double) { return thermo_property; };
    SpeciesThermoModel thermo_model;
    thermo_model.gibbs_energy     = thermo_property_fn;
    thermo_model.helmholtz_energy = thermo_property_fn;
    thermo_model.internal_energy  = thermo_property_fn;
    thermo_model.enthalpy         = thermo_property_fn;
    thermo_model.entropy          = thermo_property_fn;
    thermo_model.volume           = thermo_property_fn;
    thermo_model.heat_capacity_cp = thermo_property_fn;

    std::vector<Species> species_phase1(3);

    species_phase1[0].setName("H2O");
    species_phase1[0].setElements({"H", "O"});
    species_phase1[0].setElementAtoms({2, 1});
    species_phase1[0].setThermoModel(thermo_model);

    species_phase1[1].setName("H+");
    species_phase1[1].setElements({"H"});
    species_phase1[1].setElementAtoms({1});
    species_phase1[1].setThermoModel(thermo_model);

    species_phase1[2].setName("OH-");
    species_phase1[2].setElements({"H", "O"});
    species_phase1[2].setElementAtoms({1, 1});
    species_phase1[2].setThermoModel(thermo_model);

    std::vector<Species> species_phase2(2);

    species_phase2[0].setName("CO2(g)");
    species_phase2[0].setElements({"C", "O"});
    species_phase2[0].setElementAtoms({1, 2});
    species_phase2[0].setThermoModel(thermo_model);

    species_phase2[1].setName("H2O(g)");
    species_phase2[1].setElements({"H", "O"});
    species_phase2[1].setElementAtoms({2, 1});
    species_phase2[1].setThermoModel(thermo_model);

    std::vector<Phase> phases(2);

    phases[0].setName("Aqueous");
    phases[0].setSpecies(species_phase1);

    phases[1].setName("Gaseous");
    phases[1].setSpecies(species_phase2);

    Multiphase multiphase(phases);

    return multiphase;
}

auto test_Multiphase() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(3, multiphase.elements().size());
    ASSERT(contained("H", multiphase.elements()));
    ASSERT(contained("O", multiphase.elements()));
    ASSERT(contained("C", multiphase.elements()));
    ASSERT_EQUAL(5, multiphase.species().size());
    ASSERT_EQUAL("H2O", multiphase.species()[0].name());
    ASSERT_EQUAL("H+", multiphase.species()[1].name());
    ASSERT_EQUAL("OH-", multiphase.species()[2].name());
    ASSERT_EQUAL("CO2(g)", multiphase.species()[3].name());
    ASSERT_EQUAL("H2O(g)", multiphase.species()[4].name());
    ASSERT_EQUAL("Aqueous", multiphase.phases()[0].name());
    ASSERT_EQUAL("Gaseous", multiphase.phases()[1].name());
}

auto test_numElements() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(3, numElements(multiphase));
}

auto test_numSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(5, numSpecies(multiphase));
}

auto test_numPhases() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(2, numPhases(multiphase));
}

auto test_containsElement() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT(containsElement(multiphase, "H"));
    ASSERT(containsElement(multiphase, "O"));
    ASSERT(containsElement(multiphase, "C"));
    ASSERT(not containsElement(multiphase, "N"));
    ASSERT(not containsElement(multiphase, ""));
}

auto test_containsSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT(containsSpecies(multiphase, "H2O"));
    ASSERT(containsSpecies(multiphase, "H+"));
    ASSERT(containsSpecies(multiphase, "OH-"));
    ASSERT(containsSpecies(multiphase, "CO2(g)"));
    ASSERT(containsSpecies(multiphase, "H2O(g)"));
    ASSERT(not containsSpecies(multiphase, "NO4"));
    ASSERT(not containsSpecies(multiphase, ""));
}

auto test_containsPhase() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT(containsPhase(multiphase, "Aqueous"));
    ASSERT(containsPhase(multiphase, "Gaseous"));
    ASSERT(not containsPhase(multiphase, "Mineral"));
    ASSERT(not containsPhase(multiphase, ""));
}

auto test_indexElement() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(iC, indexElement(multiphase, "C"));
    ASSERT_EQUAL(iH, indexElement(multiphase, "H"));
    ASSERT_EQUAL(iO, indexElement(multiphase, "O"));
    ASSERT_EQUAL(numElements(multiphase), indexElement(multiphase, "N"));
    ASSERT_EQUAL(numElements(multiphase), indexElement(multiphase, ""));
}

auto test_indicesElements() -> void
{
    Multiphase multiphase = createMultiphase();
    std::vector<std::string> elements1 = {"C", "H"};
    std::vector<std::string> elements2 = {"H", "O"};
    std::vector<std::string> elements3 = {"O", "C", "H"};
    std::vector<std::string> elements4 = {"N", "C", ""};
    Indices indices1 = {iC, iH};
    Indices indices2 = {iH, iO};
    Indices indices3 = {iO, iC, iH};
    Indices indices4 = {3, iC, 3};
    ASSERT_EQUAL(indices1, indicesElements(multiphase, elements1));
    ASSERT_EQUAL(indices2, indicesElements(multiphase, elements2));
    ASSERT_EQUAL(indices3, indicesElements(multiphase, elements3));
    ASSERT_EQUAL(indices4, indicesElements(multiphase, elements4));
}

auto test_indexSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(iH2O,  indexSpecies(multiphase, "H2O"));
    ASSERT_EQUAL(iHp,   indexSpecies(multiphase, "H+"));
    ASSERT_EQUAL(iOHm,  indexSpecies(multiphase, "OH-"));
    ASSERT_EQUAL(iCO2g, indexSpecies(multiphase, "CO2(g)"));
    ASSERT_EQUAL(iH2Og, indexSpecies(multiphase, "H2O(g)"));
    ASSERT_EQUAL(numSpecies(multiphase), indexSpecies(multiphase, "NH4(g)"));
    ASSERT_EQUAL(numSpecies(multiphase), indexSpecies(multiphase, ""));
}

auto test_indicesSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    std::vector<std::string> species1 = {"CO2(g)", "H+"};
    std::vector<std::string> species2 = {"H2O", "OH-"};
    std::vector<std::string> species3 = {"CO(g)", "H2O(g)", ""};
    Indices indices1 = {iCO2g, iHp};
    Indices indices2 = {iH2O, iOHm};
    Indices indices3 = {5, iH2Og, 5};
    ASSERT_EQUAL(indices1, indicesSpecies(multiphase, species1));
    ASSERT_EQUAL(indices2, indicesSpecies(multiphase, species2));
    ASSERT_EQUAL(indices3, indicesSpecies(multiphase, species3));
}

auto test_indexPhase() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(0, indexPhase(multiphase, "Aqueous"));
    ASSERT_EQUAL(1, indexPhase(multiphase, "Gaseous"));
    ASSERT_EQUAL(numPhases(multiphase), indexPhase(multiphase, "Mineral"));
    ASSERT_EQUAL(numPhases(multiphase), indexPhase(multiphase, ""));
}

auto test_indicesPhases() -> void
{
    Multiphase multiphase = createMultiphase();
    std::vector<std::string> phases1 = {"Aqueous", "Gaseous"};
    std::vector<std::string> phases2 = {"Mineral", "", "Gaseous"};
    Indices indices1 = {0, 1};
    Indices indices2 = {2, 2, 1};
    ASSERT_EQUAL(indices1, indicesPhases(multiphase, phases1));
    ASSERT_EQUAL(indices2, indicesPhases(multiphase, phases2));
}

auto test_indexBeginSpeciesInPhase() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(0, indexBeginSpeciesInPhase(multiphase, 0));
    ASSERT_EQUAL(3, indexBeginSpeciesInPhase(multiphase, 1));
    ASSERT_EQUAL(numSpecies(multiphase), indexBeginSpeciesInPhase(multiphase, 2));
}

auto test_indexEndSpeciesInPhase() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(3, indexEndSpeciesInPhase(multiphase, 0));
    ASSERT_EQUAL(5, indexEndSpeciesInPhase(multiphase, 1));
    ASSERT_EQUAL(numSpecies(multiphase), indexEndSpeciesInPhase(multiphase, 2));
}

auto test_indicesElementsInSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices indices1 = {iH, iO};
    Indices indices2 = {iC, iO};
    ASSERT(equal(indices1, indicesElementsInSpecies(multiphase, iH2O)));
    ASSERT(equal(indices2, indicesElementsInSpecies(multiphase, iCO2g)));
}

auto test_indicesElementsInSpeciesSet() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices ispecies = {iH2O, iCO2g};
    Indices ielements = {iH, iO, iC};
    ASSERT(equal(ielements, indicesElementsInSpecies(multiphase, ispecies)));
}

auto test_indicesSpeciesInPhase() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices indices1 = {0, 1, 2};
    Indices indices2 = {3, 4};
    ASSERT(equal(indices1, indicesSpeciesInPhase(multiphase, 0)));
    ASSERT(equal(indices2, indicesSpeciesInPhase(multiphase, 1)));
}

auto test_indicesSpeciesWithElement() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices indices_with_H = {iH2O, iHp, iOHm, iH2Og};
    Indices indices_with_O = {iH2O, iOHm, iH2Og, iCO2g};
    Indices indices_with_C = {iCO2g};
    ASSERT(equal(indices_with_H, indicesSpeciesWithElement(multiphase, iH)));
    ASSERT(equal(indices_with_O, indicesSpeciesWithElement(multiphase, iO)));
    ASSERT(equal(indices_with_C, indicesSpeciesWithElement(multiphase, iC)));
}

auto test_indexPhaseWithSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(0, indexPhaseWithSpecies(multiphase, iH2O));
    ASSERT_EQUAL(0, indexPhaseWithSpecies(multiphase, iHp));
    ASSERT_EQUAL(0, indexPhaseWithSpecies(multiphase, iOHm));
    ASSERT_EQUAL(1, indexPhaseWithSpecies(multiphase, iCO2g));
    ASSERT_EQUAL(1, indexPhaseWithSpecies(multiphase, iH2Og));
}

auto test_indicesPhasesWithSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices ispecies1 = {iH2O, iOHm};
    Indices ispecies2 = {iHp, iCO2g};
    Indices ispecies3 = {iH2Og, iH2O};
    Indices iphases1 = {0};
    Indices iphases2 = {0, 1};
    Indices iphases3 = {1, 0};
    ASSERT(equal(iphases1, indicesPhasesWithSpecies(multiphase, ispecies1)));
    ASSERT(equal(iphases2, indicesPhasesWithSpecies(multiphase, ispecies2)));
    ASSERT(equal(iphases3, indicesPhasesWithSpecies(multiphase, ispecies3)));
}

auto test_localIndexSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    ASSERT_EQUAL(0, localIndexSpecies(multiphase, iH2O));
    ASSERT_EQUAL(1, localIndexSpecies(multiphase, iHp));
    ASSERT_EQUAL(2, localIndexSpecies(multiphase, iOHm));
    ASSERT_EQUAL(0, localIndexSpecies(multiphase, iCO2g));
    ASSERT_EQUAL(1, localIndexSpecies(multiphase, iH2Og));
}

auto test_indexMapSpeciesToElements() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices ielementsH2O  = {iH, iO};
    Indices ielementsHp   = {iH};
    Indices ielementsOHm  = {iH, iO};
    Indices ielementsCO2g = {iC, iO};
    Indices ielementsH2Og = {iH, iO};
    auto map = indexMapSpeciesToElements(multiphase);
    ASSERT(equal(ielementsH2O,  map[iH2O]));
    ASSERT(equal(ielementsHp,   map[iHp]));
    ASSERT(equal(ielementsOHm,  map[iOHm]));
    ASSERT(equal(ielementsCO2g, map[iCO2g]));
    ASSERT(equal(ielementsH2Og, map[iH2Og]));
}

auto test_indexMapElementToSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices ispeciesH = {iH2O, iHp, iOHm, iH2Og};
    Indices ispeciesO = {iH2O, iOHm, iH2Og, iCO2g};
    Indices ispeciesC = {iCO2g};
    auto map = indexMapElementToSpecies(multiphase);
    ASSERT(equal(ispeciesH, map[iH]));
    ASSERT(equal(ispeciesO, map[iO]));
    ASSERT(equal(ispeciesC, map[iC]));
}

auto test_indexMapPhaseToSpecies() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices ispecies_aqueous = {0, 1, 2};
    Indices ispecies_gaseous = {3, 4};
    auto map = indexMapPhaseToSpecies(multiphase);
    ASSERT_EQUAL(ispecies_aqueous, map[0]);
    ASSERT_EQUAL(ispecies_gaseous, map[1]);

}

auto test_indexMapSpeciesToPhase() -> void
{
    Multiphase multiphase = createMultiphase();
    Indices iphases = {0, 0, 0, 1, 1};
    auto map = indexMapSpeciesToPhase(multiphase);
    ASSERT_EQUAL(iphases, map);
}

#define ASSERT_EQUAL_VECTOR(estimated, actual) ASSERT(arma::all(estimated == actual))
#define ASSERT_EQUAL_VECTOR_DELTA(estimated, actual, delta) ASSERT(arma::norm(estimated - actual) < delta)
#define ASSERT_EQUAL_MATRIX(estimated, actual) ASSERT(arma::all(arma::all(estimated == actual)))
#define ASSERT_EQUAL_MATRIX_DELTA(estimated, actual, delta) ASSERT(arma::norm(estimated - actual) < delta)

auto test_formulaMatrix() -> void
{
    Multiphase multiphase = createMultiphase();
    Matrix formula_matrix(3, 5); // H2O, H+, OH-, CO2(g), H2O(g)
    formula_matrix.row(iH) = Vector{2,   1,  1,   0,      2}.t();
    formula_matrix.row(iO) = Vector{1,   0,  1,   2,      1}.t();
    formula_matrix.row(iC) = Vector{0,   0,  0,   1,      0}.t();
    Matrix formula_matrix_actual = formulaMatrix(multiphase);
    ASSERT_EQUAL_MATRIX(formula_matrix, formula_matrix_actual);
}

auto test_blockVector() -> void
{
    Multiphase multiphase = createMultiphase();
    Vector n  = {1, 2, 3, 4, 5};
    Vector n0 = {1, 2, 3};
    Vector n1 = {4, 5};
    ASSERT_EQUAL_VECTOR(n0, block(multiphase, 0, n));
    ASSERT_EQUAL_VECTOR(n1, block(multiphase, 1, n));
}

auto test_blockMatrix() -> void
{
    Multiphase multiphase = createMultiphase();
    Matrix m0 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix m1 = {1, 2, 3, 4};
    m0.reshape(3, 3);
    m1.reshape(2, 2);
    Matrix m = zeros(5, 5);
    m.submat(0, 0, 2, 2) = m0;
    m.submat(3, 3, 4, 4) = m1;
    ASSERT_EQUAL_MATRIX(m0, block(multiphase, 0, m));
    ASSERT_EQUAL_MATRIX(m1, block(multiphase, 1, m));
}

auto test_multiphaseSpeciesThermoVector() -> void
{
    Multiphase multiphase = createMultiphase();

    ThermoVector thermo_properties(1.0*ones(5), 2.0*ones(5), 3.0*ones(5));

    ASSERT_EQUAL(thermo_properties, enthalpies(multiphase, 300, 1));
    ASSERT_EQUAL(thermo_properties, entropies(multiphase, 300, 1));
    ASSERT_EQUAL(thermo_properties, gibbsEnergies(multiphase, 300, 1));
    ASSERT_EQUAL(thermo_properties, heatCapacitiesCp(multiphase, 300, 1));
    ASSERT_EQUAL(thermo_properties, helmholtzEnergies(multiphase, 300, 1));
    ASSERT_EQUAL(thermo_properties, internalEnergies(multiphase, 300, 1));
    ASSERT_EQUAL(thermo_properties, volumes(multiphase, 300, 1));
}

auto test_molarFractions() -> void
{
    Multiphase multiphase = createMultiphase();
    Vector n = {1.0, 3.0, 6.0, 3.0, 7.0};
    Vector x_val0 = {0.1, 0.3, 0.6};
    Vector x_val1 = {0.3, 0.7};

    Matrix x_ddn0 = { 0.09, -0.03, -0.06,
                     -0.01,  0.07, -0.06,
                     -0.01, -0.03,  0.04};
    Matrix x_ddn1 = { 0.07, -0.07,
                     -0.03,  0.03};
    x_ddn0.reshape(3, 3);
    x_ddn1.reshape(2, 2);

    Vector x_val(5);
    x_val.subvec(0, 2) = x_val0;
    x_val.subvec(3, 4) = x_val1;
    Vector x_ddt = zeros(5);
    Vector x_ddp = zeros(5);
    Matrix x_ddn = zeros(5, 5);
    x_ddn.submat(0, 0, 2, 2) = x_ddn0;
    x_ddn.submat(3, 3, 4, 4) = x_ddn1;

    ChemicalVector x_actual = molarFractions(multiphase, n);
    const double eps = 1.0e-16;
    ASSERT_EQUAL_VECTOR_DELTA(x_actual.val(), x_val, eps);
    ASSERT_EQUAL_VECTOR_DELTA(x_actual.ddt(), x_ddt, eps);
    ASSERT_EQUAL_VECTOR_DELTA(x_actual.ddp(), x_ddp, eps);
    ASSERT_EQUAL_MATRIX_DELTA(x_actual.ddn(), x_ddn, eps);
}

auto test_concentrations() -> void
{
    Multiphase multiphase = createMultiphase();
    // todo
}

auto test_activities() -> void
{
    Multiphase multiphase = createMultiphase();
    // todo
}

auto test_densities() -> void
{
    Multiphase multiphase = createMultiphase();
    // todo
}

} // namespace

auto testSuiteMultiphase() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_Multiphase);
    s += CUTE(test_numElements);
    s += CUTE(test_numSpecies);
    s += CUTE(test_numPhases);
    s += CUTE(test_containsElement);
    s += CUTE(test_containsSpecies);
    s += CUTE(test_containsPhase);
    s += CUTE(test_indexElement);
    s += CUTE(test_indicesElements);
    s += CUTE(test_indexSpecies);
    s += CUTE(test_indicesSpecies);
    s += CUTE(test_indexPhase);
    s += CUTE(test_indicesPhases);
    s += CUTE(test_indexBeginSpeciesInPhase);
    s += CUTE(test_indexEndSpeciesInPhase);
    s += CUTE(test_indicesElementsInSpecies);
    s += CUTE(test_indicesElementsInSpeciesSet);
    s += CUTE(test_indicesSpeciesInPhase);
    s += CUTE(test_indicesSpeciesWithElement);
    s += CUTE(test_indexPhaseWithSpecies);
    s += CUTE(test_indicesPhasesWithSpecies);
    s += CUTE(test_localIndexSpecies);
    s += CUTE(test_indexMapSpeciesToElements);
    s += CUTE(test_indexMapElementToSpecies);
    s += CUTE(test_indexMapPhaseToSpecies);
    s += CUTE(test_indexMapSpeciesToPhase);
    s += CUTE(test_formulaMatrix);
    s += CUTE(test_blockVector);
    s += CUTE(test_blockMatrix);
    s += CUTE(test_multiphaseSpeciesThermoVector);
    s += CUTE(test_molarFractions);
    s += CUTE(test_concentrations);
    s += CUTE(test_activities);
    s += CUTE(test_densities);

    return s;
}

} // namespace Reaktor
