// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "TestPartition.hpp"

// Reaktoro includes
#include <Reaktoro/Reaktoro.hpp>

namespace Reaktoro {

auto createChemicalSystem() -> ChemicalSystem
{
    ElementData element_data;
    SpeciesData species_data;

    element_data.name = "H";
    element_data.molar_mass = 1.0;
    Element eH(element_data);
 
    element_data.name = "O";
    element_data.molar_mass = 16.0;
    Element eO(element_data);

    element_data.name = "C";
    element_data.molar_mass = 12.0;
    Element eC(element_data);
 
    std::vector<Species> aqueous_species;

    species_data.name = "H2O";
    species_data.elements = {eH, eO};
    species_data.atoms = {2, 1};
    aqueous_species.emplace_back(species_data);

    species_data.name = "H+";
    species_data.elements = {eH};
    species_data.atoms = {1};
    aqueous_species.emplace_back(species_data);

    species_data.name = "OH-";
    species_data.elements = {eH, eO};
    species_data.atoms = {1, 1};
    aqueous_species.emplace_back(species_data);

    std::vector<Species> gaseous_species(2);
    gaseous_species[0].setName("H2O(g)");
    gaseous_species[0].setElements({"H", "O"});
    gaseous_species[0].setatoms({2, 1});

    gaseous_species[1].setName("CO2(g)");
    gaseous_species[1].setElements({"C", "O"});
    gaseous_species[1].setatoms({1, 2});

    std::vector<Phase> phases(2);
    phases[0].setName("Aqueous");
    phases[1].setName("Gaseous");
    phases[0].setSpecies(aqueous_species);
    phases[1].setSpecies(gaseous_species);
    ChemicalSystem system(phases);
    return system;
}

auto test_Partition() -> void
{
    Indices iequilibrium = {0, 1, 2};
    Indices ikinetic = {3, 4};
    Indices iinert = {5};
    Partition partition(iequilibrium, ikinetic, iinert);
    ASSERT_EQUAL(iequilibrium, partition.indicesEquilibriumSpecies());
    ASSERT_EQUAL(ikinetic, partition.indicesKineticSpecies());
    ASSERT_EQUAL(iinert, partition.indicesInertSpecies());
}

auto test_allEquilibrium() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition = Partition::allEquilibrium(system);
    Indices iequilibrium = range(5u);
    ASSERT_EQUAL(iequilibrium, partition.indicesEquilibriumSpecies());
    ASSERT(partition.indicesKineticSpecies().empty());
    ASSERT(partition.indicesInertSpecies().empty());
}

auto test_allKinetic() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition = Partition::allKinetic(system);
    Indices ikinetic = range(5u);
    ASSERT_EQUAL(ikinetic, partition.indicesKineticSpecies());
    ASSERT(partition.indicesEquilibriumSpecies().empty());
    ASSERT(partition.indicesInertSpecies().empty());
}

auto test_allEquilibriumExcept() -> void
{
    Indices iequilibrium = {0, 1, 3};
    Indices ikinetic = {2, 4};
    Indices iinert = {5};
    ChemicalSystem system = createChemicalSystem();
    Partition partition = Partition::allEquilibriumExcept(system, ikinetic, iinert);
    ASSERT_EQUAL(iequilibrium, partition.indicesEquilibriumSpecies());
    ASSERT_EQUAL(ikinetic, partition.indicesKineticSpecies());
    ASSERT_EQUAL(iinert, partition.indicesInertSpecies());
}

auto test_allKineticExcept() -> void
{
    Indices iequilibrium = {0, 1, 3};
    Indices ikinetic = {2, 4};
    Indices iinert = {5};
    ChemicalSystem system = createChemicalSystem();
    Partition partition = Partition::allKineticExcept(system, iequilibrium, iinert);
    ASSERT_EQUAL(iequilibrium, partition.indicesEquilibriumSpecies());
    ASSERT_EQUAL(ikinetic, partition.indicesKineticSpecies());
    ASSERT_EQUAL(iinert, partition.indicesInertSpecies());
}

auto test_numSpecies() -> void
{
    Partition partition({0,1,2}, {3,4}, {5});
    ASSERT_EQUAL(6, numSpecies(partition));
}

auto test_numEquilibriumSpecies() -> void
{
    Partition partition({0,1,2}, {3,4}, {5});
    ASSERT_EQUAL(3, numEquilibriumSpecies(partition));
}

auto test_numKineticSpecies() -> void
{
    Partition partition({0,1,2}, {3,4}, {5});
    ASSERT_EQUAL(2, numKineticSpecies(partition));
}

auto test_numInertSpecies() -> void
{
    Partition partition({0,1,2}, {3,4}, {5});
    ASSERT_EQUAL(1, numInertSpecies(partition));
}

auto test_elementIndicesInEquilibriumSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    const Index iH = indexElement(system, "H");
    const Index iO = indexElement(system, "O");
    Indices expected = {iH, iO};
    Indices actual = elementIndicesInEquilibriumSpecies(system, partition);
    ASSERT(equal(expected, actual));
}

auto test_elementIndicesInKineticSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    const Index iO = indexElement(system, "O");
    const Index iC = indexElement(system, "C");
    Indices expected = {iC, iO};
    Indices actual = elementIndicesInKineticSpecies(system, partition);
    ASSERT(equal(expected, actual));
}

auto test_elementIndicesInInertSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    const Index iH = indexElement(system, "H");
    const Index iO = indexElement(system, "O");
    Indices expected = {iH, iO};
    Indices actual = elementIndicesInInertSpecies(system, partition);
    ASSERT(equal(expected, actual));
}

auto test_phaseIndicesWithEquilibriumSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    Indices expected = {0};
    Indices actual = phaseIndicesWithEquilibriumSpecies(system, partition);
    ASSERT(equal(expected, actual));
}

auto test_phaseIndicesWithKineticSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    Indices expected = {1};
    Indices actual = phaseIndicesWithKineticSpecies(system, partition);
    ASSERT(equal(expected, actual));
}

auto test_phaseIndicesWithInertSpecies() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    Indices expected = {1};
    Indices actual = phaseIndicesWithInertSpecies(system, partition);
    ASSERT(equal(expected, actual));
}

#define ASSERT_EQUAL_ARMA(expected, actual) ASSERT(arma::norm(expected - actual) < 1.e-16)
#define ASSERT_EQUAL_MATRIX(expected, actual) ASSERT(arma::norm(expected - actual) < 1.e-16)

auto test_equilibriumRows() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Vector vec = {2.0, 3.0, 4.0, 5.0, 6.0};
    Vector expected = {2.0, 3.0, 4.0};
    Vector actual = equilibriumRows(partition, vec);
    ASSERT_EQUAL_ARMA(expected, actual);
}

auto test_kineticRows() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Vector vec = {2.0, 3.0, 4.0, 5.0, 6.0};
    Vector expected = {6.0};
    Vector actual = kineticRows(partition, vec);
    ASSERT_EQUAL_ARMA(expected, actual);
}

auto test_inertRows() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Vector vec = {2.0, 3.0, 4.0, 5.0, 6.0};
    Vector expected = {5.0};
    Vector actual = inertRows(partition, vec);
    ASSERT_EQUAL_ARMA(expected, actual);
}

auto test_equilibriumCols() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Matrix mat = arma::randn<Matrix>(3, 5);
    Matrix expected = mat.cols(arma::uvec{0,1,2});
    Matrix actual = equilibriumCols(partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto test_kineticCols() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Matrix mat = arma::randn<Matrix>(3, 5);
    Matrix expected = mat.cols(arma::uvec{4});
    Matrix actual = kineticCols(partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto test_inertCols() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Matrix mat = arma::randn<Matrix>(3, 5);
    Matrix expected = mat.cols(arma::uvec{3});
    Matrix actual = inertCols(partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto test_equilibriumRowsCols() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Matrix mat = arma::randn<Matrix>(5, 5);
    Matrix expected = mat.submat(arma::uvec{0,1,2}, arma::uvec{0,1,2});
    Matrix actual = equilibriumRowsCols(partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto test_kineticRowsCols() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Matrix mat = arma::randn<Matrix>(5, 5);
    Matrix expected = mat.submat(arma::uvec{4}, arma::uvec{4});
    Matrix actual = kineticRowsCols(partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto test_inertRowsCols() -> void
{
    Partition partition({0,1,2}, {4}, {3});
    Matrix mat = arma::randn<Matrix>(5, 5);
    Matrix expected = mat.submat(arma::uvec{3}, arma::uvec{3});
    Matrix actual = inertRowsCols(partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto test_equilibriumFormulaMatrix() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    const Index iH = indexElement(system, "H");
    const Index iO = indexElement(system, "O");
    Matrix mat = formulaMatrix(system);
    Matrix expected = mat.submat(arma::uvec{iH, iO}, arma::uvec{0,1,2});
    Matrix actual = equilibriumFormulaMatrix(system, partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto test_kineticFormulaMatrix() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    const Index iC = indexElement(system, "C");
    const Index iO = indexElement(system, "O");
    Matrix mat = formulaMatrix(system);
    Matrix expected = mat.submat(arma::uvec{iC, iO}, arma::uvec{4});
    Matrix actual = kineticFormulaMatrix(system, partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto test_inertFormulaMatrix() -> void
{
    ChemicalSystem system = createChemicalSystem();
    Partition partition({0,1,2}, {4}, {3});
    const Index iH = indexElement(system, "H");
    const Index iO = indexElement(system, "O");
    Matrix mat = formulaMatrix(system);
    Matrix expected = mat.submat(arma::uvec{iH, iO}, arma::uvec{3});
    Matrix actual = inertFormulaMatrix(system, partition, mat);
    ASSERT_EQUAL_MATRIX(expected, actual);
}

auto testSuitePartition() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_Partition);
    s += CUTE(test_allEquilibrium);
    s += CUTE(test_allKinetic);
    s += CUTE(test_allEquilibriumExcept);
    s += CUTE(test_allKineticExcept);
    s += CUTE(test_numSpecies);
    s += CUTE(test_numEquilibriumSpecies);
    s += CUTE(test_numKineticSpecies);
    s += CUTE(test_numInertSpecies);
    s += CUTE(test_elementIndicesInEquilibriumSpecies);
    s += CUTE(test_elementIndicesInKineticSpecies);
    s += CUTE(test_elementIndicesInInertSpecies);
    s += CUTE(test_phaseIndicesWithEquilibriumSpecies);
    s += CUTE(test_phaseIndicesWithKineticSpecies);
    s += CUTE(test_phaseIndicesWithInertSpecies);
    s += CUTE(test_equilibriumRows);
    s += CUTE(test_kineticRows);
    s += CUTE(test_inertRows);
    s += CUTE(test_equilibriumCols);
    s += CUTE(test_kineticCols);
    s += CUTE(test_inertCols);
    s += CUTE(test_equilibriumRowsCols);
    s += CUTE(test_kineticRowsCols);
    s += CUTE(test_inertRowsCols);
    s += CUTE(test_equilibriumFormulaMatrix);
    s += CUTE(test_kineticFormulaMatrix);
    s += CUTE(test_inertFormulaMatrix);

    return s;
}

} // namespace Reaktoro
