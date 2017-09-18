// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

namespace Reaktoro {
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

auto speciesMoles() -> Vector
{
    Vector res(5); res << 1.0, 3.0, 6.0, 3.0, 7.0;
    return res;
}

auto aqueousSpeciesMolarFractions() -> ChemicalVector
{
    Vector x_val(3); x_val << 0.1, 0.3, 0.6;
    Vector x_ddt = zeros(3);
    Vector x_ddp = zeros(3);
    Matrix x_ddn(3, 3); x_ddn << 0.09, -0.03, -0.06,
                                -0.01,  0.07, -0.06,
                                -0.01, -0.03,  0.04;
    return ChemicalVector{x_val, x_ddt, x_ddp, x_ddn};
}

auto gaseousSpeciesMolarFractions() -> ChemicalVector
{
    Vector x_val(2); x_val << 0.3, 0.7;
    Vector x_ddt = zeros(2);
    Vector x_ddp = zeros(2);
    Matrix x_ddn(2, 2); x_ddn << 0.07, -0.07,
                                -0.03,  0.03;
    return ChemicalVector{x_val, x_ddt, x_ddp, x_ddn};
}

auto aqueousPhaseDensity() -> ChemicalScalar
{
    return ChemicalScalar{1000.0, -100.0, +200.0, Vector{0.2, 0.3, 0.5}};
}

auto gaseousPhaseDensity() -> ChemicalScalar
{
    return ChemicalScalar{10.0, -1000.0, +20.0, Vector{0.6, 0.2}};
}

auto aqueousPhaseThermoModel() -> PhaseThermoModel
{
    PhaseThermoModel phase_thermo_model;
    phase_thermo_model.concentration = [=](const Vector&) { return aqueousSpeciesMolarFractions(); };
    phase_thermo_model.activity = [=](Temperature, Pressure, const Vector&) { return aqueousSpeciesMolarFractions(); };
    phase_thermo_model.density = [=](Temperature, Pressure, const Vector&) { return aqueousPhaseDensity(); };
    return phase_thermo_model;
}

auto gaseousPhaseThermoModel() -> PhaseThermoModel
{
    PhaseThermoModel phase_thermo_model;
    phase_thermo_model.concentration = [=](const Vector&) { return gaseousSpeciesMolarFractions(); };
    phase_thermo_model.activity = [=](Temperature, Pressure, const Vector&) { return gaseousSpeciesMolarFractions(); };
    phase_thermo_model.density = [=](Temperature, Pressure, const Vector&) { return gaseousPhaseDensity(); };
    return phase_thermo_model;
}

auto speciesMolarFractions() -> ChemicalVector
{
    Vector x_val(5); x_val << 0.1, 0.3, 0.6, 0.3, 0.7;
    Vector x_ddt = zeros(5);
    Vector x_ddp = zeros(5);
    Matrix x_ddn(5, 5); x_ddn << 0.09, -0.03, -0.06,  0.00,  0.00,
                                -0.01,  0.07, -0.06,  0.00,  0.00,
                                -0.01, -0.03,  0.04,  0.00,  0.00,
                                 0.00,  0.00,  0.00,  0.07, -0.07,
                                 0.00,  0.00,  0.00, -0.03,  0.03;
    return ChemicalVector{x_val, x_ddt, x_ddp, x_ddn};
}

auto phaseDensities() -> ChemicalVector
{
    ChemicalVector d(2, 5);
    d.row(0) = ChemicalScalar{1000.0, -100.0, +200.0, Vector{0.2, 0.3, 0.5, 0.0, 0.0}};
    d.row(1) = ChemicalScalar{10.0, -1000.0, +20.0, Vector{0.0, 0.0, 0.0, 0.6, 0.2}};
    return d;
}

auto createChemicalSystem() -> ChemicalSystem
{
    ThermoScalar thermo_scalar(1.0, 2.0, 3.0);
    ThermoScalarFunction thermo_scalar_fn = [=](double,double) { return thermo_scalar; };
    SpeciesThermoModel species_thermo_model;
    species_thermo_model.gibbs_energy     = thermo_scalar_fn;
    species_thermo_model.helmholtz_energy = thermo_scalar_fn;
    species_thermo_model.internal_energy  = thermo_scalar_fn;
    species_thermo_model.enthalpy         = thermo_scalar_fn;
    species_thermo_model.entropy          = thermo_scalar_fn;
    species_thermo_model.volume           = thermo_scalar_fn;
    species_thermo_model.heat_capacity = thermo_scalar_fn;

    std::vector<Species> aqueous_species(3);

    aqueous_species[0].setName("H2O");
    aqueous_species[0].setElements({"H", "O"});
    aqueous_species[0].setElementAtoms({2, 1});
    aqueous_species[0].setThermoModel(species_thermo_model);

    aqueous_species[1].setName("H+");
    aqueous_species[1].setElements({"H"});
    aqueous_species[1].setElementAtoms({1});
    aqueous_species[1].setThermoModel(species_thermo_model);

    aqueous_species[2].setName("OH-");
    aqueous_species[2].setElements({"H", "O"});
    aqueous_species[2].setElementAtoms({1, 1});
    aqueous_species[2].setThermoModel(species_thermo_model);

    std::vector<Species> gaseous_species(2);

    gaseous_species[0].setName("CO2(g)");
    gaseous_species[0].setElements({"C", "O"});
    gaseous_species[0].setElementAtoms({1, 2});
    gaseous_species[0].setThermoModel(species_thermo_model);

    gaseous_species[1].setName("H2O(g)");
    gaseous_species[1].setElements({"H", "O"});
    gaseous_species[1].setElementAtoms({2, 1});
    gaseous_species[1].setThermoModel(species_thermo_model);

    std::vector<Phase> phases(2);

    phases[0].setName("Aqueous");
    phases[0].setSpecies(aqueous_species);
    phases[0].setThermoModel(aqueousPhaseThermoModel());

    phases[1].setName("Gaseous");
    phases[1].setSpecies(gaseous_species);
    phases[1].setThermoModel(gaseousPhaseThermoModel());

    ChemicalSystem system(phases);

    return system;
}

auto test_ChemicalSystem() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(3, system.elements().size());
    ASSERT(contained("H", system.elements()));
    ASSERT(contained("O", system.elements()));
    ASSERT(contained("C", system.elements()));
    ASSERT_EQUAL(5, system.species().size());
    ASSERT_EQUAL("H2O", system.species()[0].name());
    ASSERT_EQUAL("H+", system.species()[1].name());
    ASSERT_EQUAL("OH-", system.species()[2].name());
    ASSERT_EQUAL("CO2(g)", system.species()[3].name());
    ASSERT_EQUAL("H2O(g)", system.species()[4].name());
    ASSERT_EQUAL("Aqueous", system.phases()[0].name());
    ASSERT_EQUAL("Gaseous", system.phases()[1].name());
}

auto test_numElements() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(3, numElements(system));
}

auto test_numSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(5, numSpecies(system));
}

auto test_numPhases() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(2, numPhases(system));
}

auto test_containsElement() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT(containsElement(system, "H"));
    ASSERT(containsElement(system, "O"));
    ASSERT(containsElement(system, "C"));
    ASSERT(not containsElement(system, "N"));
    ASSERT(not containsElement(system, ""));
}

auto test_containsSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT(containsSpecies(system, "H2O"));
    ASSERT(containsSpecies(system, "H+"));
    ASSERT(containsSpecies(system, "OH-"));
    ASSERT(containsSpecies(system, "CO2(g)"));
    ASSERT(containsSpecies(system, "H2O(g)"));
    ASSERT(not containsSpecies(system, "NO4"));
    ASSERT(not containsSpecies(system, ""));
}

auto test_containsPhase() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT(containsPhase(system, "Aqueous"));
    ASSERT(containsPhase(system, "Gaseous"));
    ASSERT(not containsPhase(system, "Mineral"));
    ASSERT(not containsPhase(system, ""));
}

auto test_indexElement() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(iC, indexElement(system, "C"));
    ASSERT_EQUAL(iH, indexElement(system, "H"));
    ASSERT_EQUAL(iO, indexElement(system, "O"));
    ASSERT_EQUAL(numElements(system), indexElement(system, "N"));
    ASSERT_EQUAL(numElements(system), indexElement(system, ""));
}

auto test_elementIndices() -> void
{
    ChemicalSystem system = createChemicalSystem();
    std::vector<std::string> elements1 = {"C", "H"};
    std::vector<std::string> elements2 = {"H", "O"};
    std::vector<std::string> elements3 = {"O", "C", "H"};
    std::vector<std::string> elements4 = {"N", "C", ""};
    Indices indices1 = {iC, iH};
    Indices indices2 = {iH, iO};
    Indices indices3 = {iO, iC, iH};
    Indices indices4 = {3, iC, 3};
    ASSERT_EQUAL(indices1, elementIndices(system, elements1));
    ASSERT_EQUAL(indices2, elementIndices(system, elements2));
    ASSERT_EQUAL(indices3, elementIndices(system, elements3));
    ASSERT_EQUAL(indices4, elementIndices(system, elements4));
}

auto test_elementIndicesInSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices indices1 = {iH, iO};
    Indices indices2 = {iC, iO};
    ASSERT(equal(indices1, elementIndicesInSpecies(system, iH2O)));
    ASSERT(equal(indices2, elementIndicesInSpecies(system, iCO2g)));
}

auto test_elementIndicesInSpeciesArray() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices ispecies = {iH2O, iCO2g};
    Indices ielements = {iH, iO, iC};
    ASSERT(equal(ielements, elementIndicesInSpecies(system, ispecies)));
}

auto test_indexSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(iH2O,  indexSpecies(system, "H2O"));
    ASSERT_EQUAL(iHp,   indexSpecies(system, "H+"));
    ASSERT_EQUAL(iOHm,  indexSpecies(system, "OH-"));
    ASSERT_EQUAL(iCO2g, indexSpecies(system, "CO2(g)"));
    ASSERT_EQUAL(iH2Og, indexSpecies(system, "H2O(g)"));
    ASSERT_EQUAL(numSpecies(system), indexSpecies(system, "NH4(g)"));
    ASSERT_EQUAL(numSpecies(system), indexSpecies(system, ""));
}

auto test_speciesIndices() -> void
{
    ChemicalSystem system = createChemicalSystem();
    std::vector<std::string> species1 = {"CO2(g)", "H+"};
    std::vector<std::string> species2 = {"H2O", "OH-"};
    std::vector<std::string> species3 = {"CO(g)", "H2O(g)", ""};
    Indices indices1 = {iCO2g, iHp};
    Indices indices2 = {iH2O, iOHm};
    Indices indices3 = {5, iH2Og, 5};
    ASSERT_EQUAL(indices1, speciesIndices(system, species1));
    ASSERT_EQUAL(indices2, speciesIndices(system, species2));
    ASSERT_EQUAL(indices3, speciesIndices(system, species3));
}

auto test_speciesBeginIndexInPhase() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(0, speciesBeginIndexInPhase(system, 0));
    ASSERT_EQUAL(3, speciesBeginIndexInPhase(system, 1));
    ASSERT_EQUAL(numSpecies(system), speciesBeginIndexInPhase(system, 2));
}

auto test_speciesEndIndexInPhase() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(3, speciesEndIndexInPhase(system, 0));
    ASSERT_EQUAL(5, speciesEndIndexInPhase(system, 1));
    ASSERT_EQUAL(numSpecies(system), speciesEndIndexInPhase(system, 2));
}

auto test_speciesIndicesInPhase() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices indices1 = {0, 1, 2};
    Indices indices2 = {3, 4};
    ASSERT(equal(indices1, speciesIndicesInPhase(system, 0)));
    ASSERT(equal(indices2, speciesIndicesInPhase(system, 1)));
}

auto test_speciesIndicesWithElement() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices indices_with_H = {iH2O, iHp, iOHm, iH2Og};
    Indices indices_with_O = {iH2O, iOHm, iH2Og, iCO2g};
    Indices indices_with_C = {iCO2g};
    ASSERT(equal(indices_with_H, speciesIndicesWithElement(system, iH)));
    ASSERT(equal(indices_with_O, speciesIndicesWithElement(system, iO)));
    ASSERT(equal(indices_with_C, speciesIndicesWithElement(system, iC)));
}

auto test_speciesLocalIndex() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(0, speciesLocalIndex(system, iH2O));
    ASSERT_EQUAL(1, speciesLocalIndex(system, iHp));
    ASSERT_EQUAL(2, speciesLocalIndex(system, iOHm));
    ASSERT_EQUAL(0, speciesLocalIndex(system, iCO2g));
    ASSERT_EQUAL(1, speciesLocalIndex(system, iH2Og));
}

auto test_indexPhase() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(0, indexPhase(system, "Aqueous"));
    ASSERT_EQUAL(1, indexPhase(system, "Gaseous"));
    ASSERT_EQUAL(numPhases(system), indexPhase(system, "Mineral"));
    ASSERT_EQUAL(numPhases(system), indexPhase(system, ""));
}

auto test_phaseIndices() -> void
{
    ChemicalSystem system = createChemicalSystem();
    std::vector<std::string> phases1 = {"Aqueous", "Gaseous"};
    std::vector<std::string> phases2 = {"Mineral", "", "Gaseous"};
    Indices indices1 = {0, 1};
    Indices indices2 = {2, 2, 1};
    ASSERT_EQUAL(indices1, phaseIndices(system, phases1));
    ASSERT_EQUAL(indices2, phaseIndices(system, phases2));
}

auto test_indexPhaseWithSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    ASSERT_EQUAL(0, indexPhaseWithSpecies(system, iH2O));
    ASSERT_EQUAL(0, indexPhaseWithSpecies(system, iHp));
    ASSERT_EQUAL(0, indexPhaseWithSpecies(system, iOHm));
    ASSERT_EQUAL(1, indexPhaseWithSpecies(system, iCO2g));
    ASSERT_EQUAL(1, indexPhaseWithSpecies(system, iH2Og));
}

auto test_phaseIndicesWithSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices ispecies1 = {iH2O, iOHm};
    Indices ispecies2 = {iHp, iCO2g};
    Indices ispecies3 = {iH2Og, iH2O};
    Indices iphases1 = {0};
    Indices iphases2 = {0, 1};
    Indices iphases3 = {1, 0};
    ASSERT(equal(iphases1, phaseIndicesWithSpecies(system, ispecies1)));
    ASSERT(equal(iphases2, phaseIndicesWithSpecies(system, ispecies2)));
    ASSERT(equal(iphases3, phaseIndicesWithSpecies(system, ispecies3)));
}

auto test_indexMapSpeciesToElements() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices ielementsH2O  = {iH, iO};
    Indices ielementsHp   = {iH};
    Indices ielementsOHm  = {iH, iO};
    Indices ielementsCO2g = {iC, iO};
    Indices ielementsH2Og = {iH, iO};
    auto map = indexMapSpeciesToElements(system);
    ASSERT(equal(ielementsH2O,  map[iH2O]));
    ASSERT(equal(ielementsHp,   map[iHp]));
    ASSERT(equal(ielementsOHm,  map[iOHm]));
    ASSERT(equal(ielementsCO2g, map[iCO2g]));
    ASSERT(equal(ielementsH2Og, map[iH2Og]));
}

auto test_indexMapElementToSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices ispeciesH = {iH2O, iHp, iOHm, iH2Og};
    Indices ispeciesO = {iH2O, iOHm, iH2Og, iCO2g};
    Indices ispeciesC = {iCO2g};
    auto map = indexMapElementToSpecies(system);
    ASSERT(equal(ispeciesH, map[iH]));
    ASSERT(equal(ispeciesO, map[iO]));
    ASSERT(equal(ispeciesC, map[iC]));
}

auto test_indexMapPhaseToSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices ispecies_aqueous = {0, 1, 2};
    Indices ispecies_gaseous = {3, 4};
    auto map = indexMapPhaseToSpecies(system);
    ASSERT_EQUAL(ispecies_aqueous, map[0]);
    ASSERT_EQUAL(ispecies_gaseous, map[1]);

}

auto test_indexMapSpeciesToPhase() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Indices iphases = {0, 0, 0, 1, 1};
    auto map = indexMapSpeciesToPhase(system);
    ASSERT_EQUAL(iphases, map);
}

#define ASSERT_EQUAL_ARMA(estimated, actual) ASSERT(arma::all(estimated == actual))
#define ASSERT_EQUAL_VECTOR_DELTA(estimated, actual, delta) ASSERT(arma::norm(estimated - actual) < delta)
#define ASSERT_EQUAL_MATRIX(estimated, actual) ASSERT(arma::all(arma::all(estimated == actual)))
#define ASSERT_EQUAL_MATRIX_DELTA(estimated, actual, delta) ASSERT(arma::norm(estimated - actual) < delta)

auto test_formulaMatrix() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Matrix formula_matrix(3, 5); // H2O, H+, OH-, CO2(g), H2O(g)
    formula_matrix.row(iH) = Vector{2,   1,  1,   0,      2}.t();
    formula_matrix.row(iO) = Vector{1,   0,  1,   2,      1}.t();
    formula_matrix.row(iC) = Vector{0,   0,  0,   1,      0}.t();
    Matrix formula_matrix_actual = formulaMatrix(system);
    ASSERT_EQUAL_MATRIX(formula_matrix, formula_matrix_actual);
}

auto test_blockVector() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Vector n  = {1, 2, 3, 4, 5};
    Vector n0 = {1, 2, 3};
    Vector n1 = {4, 5};
    ASSERT_EQUAL_ARMA(n0, block(system, 0, n));
    ASSERT_EQUAL_ARMA(n1, block(system, 1, n));
}

auto test_blockMatrix() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Matrix m0 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Matrix m1 = {1, 2, 3, 4};
    m0.reshape(3, 3);
    m1.reshape(2, 2);
    Matrix m = zeros(5, 5);
    m.submat(0, 0, 2, 2) = m0;
    m.submat(3, 3, 4, 4) = m1;
    ASSERT_EQUAL_MATRIX(m0, block(system, 0, m));
    ASSERT_EQUAL_MATRIX(m1, block(system, 1, m));
}

auto test_systemSpeciesInterpolatedThermoProperties() -> void
{
    ChemicalSystem system = createChemicalSystem();

    ThermoVector thermo_properties(1.0*ones(5), 2.0*ones(5), 3.0*ones(5));

    ASSERT_EQUAL(thermo_properties, enthalpies(system, 300, 1));
    ASSERT_EQUAL(thermo_properties, entropies(system, 300, 1));
    ASSERT_EQUAL(thermo_properties, standardGibbsEnergies(system, 300, 1));
    ASSERT_EQUAL(thermo_properties, heatCapacitiesCp(system, 300, 1));
    ASSERT_EQUAL(thermo_properties, helmholtzEnergies(system, 300, 1));
    ASSERT_EQUAL(thermo_properties, internalEnergies(system, 300, 1));
    ASSERT_EQUAL(thermo_properties, volumes(system, 300, 1));
}

auto test_moleFractions() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Vector n = speciesMoles();
    ChemicalVector x = speciesMolarFractions();
    ChemicalVector x_actual = moleFractions(system, n);
    const double eps = 1.0e-16;
    ASSERT_EQUAL_VECTOR_DELTA(x_actual.val, x.val, eps);
    ASSERT_EQUAL_VECTOR_DELTA(x_actual.ddT(), x.ddT(), eps);
    ASSERT_EQUAL_VECTOR_DELTA(x_actual.ddP(), x.ddP(), eps);
    ASSERT_EQUAL_MATRIX_DELTA(x_actual.ddn, x.ddn, eps);
}

auto test_phasesThermoModels() -> void
{
    Vector n = {2.0, 8.0};
    ChemicalScalar rho(1000.0, -100.0, +200.0, Vector{0.2, 0.3});
    ChemicalVector c(2, 2);
    c.row(0) = ChemicalScalar(0.2, 0.0, 0.0, Vector{+0.08, -0.02});
    c.row(1) = ChemicalScalar(0.8, 0.0, 0.0, Vector{-0.08, +0.02});
    PhaseThermoModel thermo_model;
    thermo_model.concentration = [=](const Vector&) { return c; };
    thermo_model.activity = [=](Temperature, Pressure, const Vector&) { return c; };
    thermo_model.density = [=](Temperature, Pressure, const Vector&) { return rho; };
    Phase phase;
    phase.setSpecies(std::vector<Species>(2));
    phase.setThermoModel(thermo_model);
    ASSERT_EQUAL(c, concentrations(phase, n));
    ASSERT_EQUAL(c, activities(phase, 300.0, 1.0, n));
    ASSERT_EQUAL(rho, density(phase, 300.0, 1.0, n));
}

auto test_concentrations() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Vector n = speciesMoles();
    ChemicalVector c = speciesMolarFractions();
    ASSERT_EQUAL(c, concentrations(system, n));
}

auto test_activities() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Vector n = speciesMoles();
    ChemicalVector a = speciesMolarFractions();
    ASSERT_EQUAL(a, activities(system, 300.0, 1.0, n));
}

auto test_densities() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Vector n = speciesMoles();
    ChemicalVector d = phaseDensities();
    ChemicalVector d_actual = densities(system, 300.0, 1.0, n);
    ASSERT_EQUAL(d, d_actual);
}

} // namespace

auto testSuiteChemicalSystem() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_ChemicalSystem);
    s += CUTE(test_numElements);
    s += CUTE(test_numSpecies);
    s += CUTE(test_numPhases);
    s += CUTE(test_containsElement);
    s += CUTE(test_containsSpecies);
    s += CUTE(test_containsPhase);
    s += CUTE(test_indexElement);
    s += CUTE(test_elementIndices);
    s += CUTE(test_elementIndicesInSpecies);
    s += CUTE(test_elementIndicesInSpeciesArray);
    s += CUTE(test_indexSpecies);
    s += CUTE(test_speciesIndices);
    s += CUTE(test_speciesBeginIndexInPhase);
    s += CUTE(test_speciesEndIndexInPhase);
    s += CUTE(test_speciesIndicesInPhase);
    s += CUTE(test_speciesIndicesWithElement);
    s += CUTE(test_speciesLocalIndex);
    s += CUTE(test_indexPhase);
    s += CUTE(test_phaseIndices);
    s += CUTE(test_indexPhaseWithSpecies);
    s += CUTE(test_phaseIndicesWithSpecies);
    s += CUTE(test_indexMapSpeciesToElements);
    s += CUTE(test_indexMapElementToSpecies);
    s += CUTE(test_indexMapPhaseToSpecies);
    s += CUTE(test_indexMapSpeciesToPhase);
    s += CUTE(test_formulaMatrix);
    s += CUTE(test_blockVector);
    s += CUTE(test_blockMatrix);
    s += CUTE(test_systemSpeciesInterpolatedThermoProperties);
    s += CUTE(test_molarFractions);
    s += CUTE(test_phasesThermoModels);
    s += CUTE(test_concentrations);
    s += CUTE(test_activities);
    s += CUTE(test_densities);

    return s;
}

} // namespace Reaktoro
