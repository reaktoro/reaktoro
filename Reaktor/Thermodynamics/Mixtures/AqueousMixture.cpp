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

#include "AqueousMixture.hpp"

// C++ includes
#include <algorithm>

// Reaktor includes
#include <Reaktor/Common/SetUtils.hpp>

namespace Reaktor {
namespace internal {

auto speciesIndex(const std::vector<AqueousSpecies>& solution, const std::string& name) -> Index
{
    for(Index i = 0; i < solution.size(); ++i)
        if(solution[i].name == name) return i;
    return solution.size();
}

auto waterIndex(const std::vector<AqueousSpecies>& solution) -> Index
{
    auto comparer = [](const AqueousSpecies& species)
    {
        return species.name == "H2O(l)";
    };

    return std::find_if(solution.begin(), solution.end(), comparer) - solution.begin();
}

auto indicesChargedSpecies(const std::vector<AqueousSpecies>& solution) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < solution.size(); ++i)
        if(solution[i].charge != 0)
            indices.push_back(i);
    return indices;
}

auto indicesNeutralSpecies(const std::vector<AqueousSpecies>& solution) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < solution.size(); ++i)
        if(solution[i].charge == 0)
            indices.push_back(i);
    return indices;
}

auto indicesCations(const std::vector<AqueousSpecies>& solution) -> Indices
{
    Indices indices_cations;
    for(unsigned i = 0; i < solution.size(); ++i)
        if(solution[i].charge > 0)
            indices_cations.push_back(i);
    return indices_cations;
}

auto indicesAnions(const std::vector<AqueousSpecies>& solution) -> Indices
{
    Indices indices_anions;
    for(unsigned i = 0; i < solution.size(); ++i)
        if(solution[i].charge < 0)
            indices_anions.push_back(i);
    return indices_anions;
}

auto chargedSpeciesLocalIndex(const std::vector<AqueousSpecies>& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return index(idx, indicesChargedSpecies(solution));
}

auto neutralSpeciesLocalIndex(const std::vector<AqueousSpecies>& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return index(idx, indicesNeutralSpecies(solution));
}

auto cationLocalIndex(const std::vector<AqueousSpecies>& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return index(idx, indicesCations(solution));
}

auto anionLocalIndex(const std::vector<AqueousSpecies>& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return index(idx, indicesAnions(solution));
}

auto speciesNames(const std::vector<AqueousSpecies>& solution) -> std::vector<std::string>
{
    const unsigned nspecies = solution.size();
    std::vector<std::string> names(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        names[i] = solution[i].name;
    return names;
}

auto chargedSpeciesNames(const std::vector<AqueousSpecies>& solution) -> std::vector<std::string>
{
    return extract(speciesNames(solution), indicesChargedSpecies(solution));
}

auto neutralSpeciesNames(const std::vector<AqueousSpecies>& solution) -> std::vector<std::string>
{
    return extract(speciesNames(solution), indicesNeutralSpecies(solution));
}

auto cationNames(const std::vector<AqueousSpecies>& solution) -> std::vector<std::string>
{
    return extract(speciesNames(solution), indicesCations(solution));
}

auto anionNames(const std::vector<AqueousSpecies>& solution) -> std::vector<std::string>
{
    return extract(speciesNames(solution), indicesAnions(solution));
}

auto chargedSpeciesCharges(const std::vector<AqueousSpecies>& solution) -> Vector
{
    return rows(speciesCharges(solution), indicesChargedSpecies(solution));
}

auto cationCharges(const std::vector<AqueousSpecies>& solution) -> Vector
{
    return rows(speciesCharges(solution), indicesCations(solution));
}

auto anionCharges(const std::vector<AqueousSpecies>& solution) -> Vector
{
    return rows(speciesCharges(solution), indicesAnions(solution));
}

auto dissociationMatrix(const std::vector<AqueousSpecies>& solution) -> Matrix
{
    // The indices of the neutral and charged species
    const Indices indices_neutral = indicesNeutralSpecies(solution);
    const Indices indices_charged = indicesChargedSpecies(solution);

    // Gets the stoichiometry of the i-th charged species in the j-th neutral species
    auto stoichiometry = [&](unsigned i, unsigned j) -> double
    {
        const Index ineutral = indices_neutral[i];
        const Index icharged = indices_charged[j];
        const AqueousSpecies& neutral = solution[ineutral];
        const AqueousSpecies& charged = solution[icharged];
        const auto iter = neutral.dissociation.find(charged.name);
        return iter != neutral.dissociation.end() ? iter->second : 0.0;
    };

    // Assemble the dissociation matrix of the neutral species with respect to the charged species
    const unsigned num_charged_species = indices_charged.size();
    const unsigned num_neutral_species = indices_neutral.size();
    Matrix dissociation_matrix = zeros(num_neutral_species, num_charged_species);
    for(unsigned i = 0; i < num_neutral_species; ++i)
        for(unsigned j = 0; j < num_charged_species; ++j)
            dissociation_matrix(i, j) = stoichiometry(i, j);
    return dissociation_matrix;
}

} // namespace internal

AqueousMixture::AqueousMixture()
{}

AqueousMixture::AqueousMixture(const std::vector<AqueousSpecies>& species)
: GeneralMixture<AqueousSpecies>(species)
{
    // Initialize the index of the water species
    idx_water = speciesIndex("H2O(l)");

    // Initialize the indices of the neutral aqueous species
    idx_neutral_species = internal::indicesNeutralSpecies(species);

    // Initialize the indices of the charged aqueous species
    idx_charged_species = internal::indicesChargedSpecies(species);

    // Initialize the indices of the cations
    idx_cations = internal::indicesCations(species);

    // Initialize the indices of the anions
    idx_anions = internal::indicesAnions(species);

    // Initialize the dissociation matrix of the neutral species w.r.t. the charged species
    dissociation_matrix = internal::dissociationMatrix(species);
}

AqueousMixture::~AqueousMixture()
{}

auto AqueousMixture::numNeutralSpecies() const -> unsigned
{
    return idx_neutral_species.size();
}

auto AqueousMixture::numChargedSpecies() const -> unsigned
{
    return idx_charged_species.size();
}

auto AqueousMixture::indicesNeutralSpecies() const -> const Indices&
{
    return idx_neutral_species;
}

auto AqueousMixture::indicesChargedSpecies() const -> const Indices&
{
    return idx_charged_species;
}

auto AqueousMixture::indicesCations() const -> const Indices&
{
    return idx_cations;
}

auto AqueousMixture::indicesAnions() const -> const Indices&
{
    return idx_anions;
}

auto AqueousMixture::indexWater() const -> Index
{
    return idx_water;
}

auto AqueousMixture::dissociationMatrix() const -> const Matrix&
{
    return dissociation_matrix;
}

auto AqueousMixture::indexNeutralSpecies(const std::string& name) const -> Index
{
    const Index idx = speciesIndex(name);
    return index(idx, idx_neutral_species);
}

auto AqueousMixture::indexChargedSpecies(const std::string& name) const -> Index
{
    const Index idx = speciesIndex(name);
    return index(idx, idx_charged_species);
}

auto AqueousMixture::indexCation(const std::string& name) const -> Index
{
    const Index idx = speciesIndex(name);
    return index(idx, idx_cations);
}

auto AqueousMixture::indexAnion(const std::string& name) const -> Index
{
    const Index idx = speciesIndex(name);
    return index(idx, idx_anions);
}

auto AqueousMixture::namesNeutralSpecies() const -> std::vector<std::string>
{
    return extract(speciesNames(), indicesNeutralSpecies());
}

auto AqueousMixture::namesChargedSpecies() const -> std::vector<std::string>
{
    return extract(speciesNames(), indicesChargedSpecies());
}

auto AqueousMixture::namesCations() const -> std::vector<std::string>
{
    return extract(speciesNames(), indicesCations());
}

auto AqueousMixture::namesAnions() const -> std::vector<std::string>
{
    return extract(speciesNames(), indicesAnions());
}

auto AqueousMixture::molalities(const Vector& n) const -> ChemicalVector
{
    const unsigned size = numSpecies();

    const double nw = n[idx_water];

    const Vector m = 55.508 * n/nw;

    Matrix dmdn = zeros(size, size);

    for(unsigned i = 0; i < size; ++i)
        dmdn(i, i) = m[i]/n[i];

    for(unsigned i = 0; i < size; ++i)
        dmdn(i, idx_water) -= m[i]/nw;

    return partialVector(m, dmdn);
}

auto AqueousMixture::stoichiometricMolalities(const ChemicalVector& m) const -> ChemicalVector
{
    // The molalities of the ionic species
    ChemicalVector m_ions;
    func(m_ions) = rows(idx_ions, func(m));
    grad(m_ions) = rows(idx_ions, grad(m));

    // The molalities of the complex species
    ChemicalVector m_complexes;
    func(m_complexes) = rows(idx_complexes, func(m));
    grad(m_complexes) = rows(idx_complexes, grad(m));

    // The stoichiometric molalities of the ionic species
    ChemicalVector ms;
    func(ms) = func(m_ions) + dissociation_matrix.transpose() * func(m_complexes);
    grad(ms) = grad(m_ions) + dissociation_matrix.transpose() * grad(m_complexes);

    return ms;
}

auto AqueousMixture::effectiveIonicStrength(const ChemicalVector& m) const -> ChemicalScalar
{
    const unsigned num_species = numSpecies();

    ChemicalScalar Ie = partialScalar(0.0, Vector(num_species));

    func(Ie) = 0.5 * (z.array() * z.array() * func(m).array()).sum();

    for(unsigned j = 0; j < num_species; ++j)
        grad(Ie)[j] = 0.5 * (z.array() * z.array() * grad(m).col(j).array()).sum();

    return Ie;
}

auto AqueousMixture::stoichiometricIonicStrength(const ChemicalVector& ms) const -> ChemicalScalar
{
    const unsigned num_species = numSpecies();

    ChemicalScalar Is = partialScalar(0.0, Vector(num_species));

    func(Is) = 0.5 * (zi.array() * zi.array() * func(ms).array()).sum();

    for(unsigned j = 0; j < num_species; ++j)
        grad(Is)[j] = 0.5 * (zi.array() * zi.array() * grad(ms).col(j).array()).sum();

    return Is;
}

} // namespace Reaktor
