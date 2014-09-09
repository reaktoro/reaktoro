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
#include <Reaktor/Common/MatrixUtils.hpp>
#include <Reaktor/Common/SetUtils.hpp>

namespace Reaktor {

AqueousMixture::AqueousMixture()
{}

AqueousMixture::AqueousMixture(const std::vector<AqueousSpecies>& species)
: GeneralMixture<AqueousSpecies>(species)
{
    // Initialize the index of the water species
    index_water = indexSpecies(*this, "H2O(l)");

    // Initialize the indices of the charged and neutral aqueous species
    for(unsigned i = 0; i < species.size(); ++i)
        if(species[i].charge != 0) indices_charged_species.push_back(i);
        else indices_neutral_species.push_back(i);

    // Initialize the electrical charges of the aqueous species
    z.resize(species.size());
    for(unsigned i = 0; i < species.size(); ++i)
        z[i] = species[i].charge;

    // Initialize the electrical charges of the ions
    z_charged = rows(indices_charged_species, z);

    // Auxiliary function that gets the stoichiometry of the i-th charged species in the j-th neutral species
    auto stoichiometry = [&](unsigned i, unsigned j) -> double
    {
    	const auto& complex = species[indices_neutral_species[i]];
    	const auto& ion = species[indices_charged_species[j]].name;
    	const Index k = find(ion, complex.dissociation.ions);
    	return k < complex.dissociation.ions.size() ?
			complex.dissociation.stoichiometries[k] : 0.0;
    };

    // Assemble the dissociation matrix of the neutral species with respect to the charged species
    const unsigned num_charged_species = indices_charged_species.size();
	const unsigned num_neutral_species = indices_neutral_species.size();
    dissociation_matrix = zeros(num_charged_species, num_neutral_species);
    for(unsigned i = 0; i < num_neutral_species; ++i)
        for(unsigned j = 0; j < num_charged_species; ++j)
            dissociation_matrix(i, j) = stoichiometry(i, j);
}

AqueousMixture::~AqueousMixture()
{}

auto AqueousMixture::charges() const -> Vector
{
    return z;
}

auto AqueousMixture::indicesChargedSpecies() const -> const Indices&
{
    return indices_charged_species;
}

auto AqueousMixture::indicesNeutralSpecies() const -> const Indices&
{
    return indices_neutral_species;
}

auto AqueousMixture::dissociationMatrix() const -> const Matrix&
{
    return dissociation_matrix;
}

auto numChargedSpecies(const AqueousMixture& mixture) -> unsigned
{
	return mixture.indicesChargedSpecies().size();
}

auto numNeutralSpecies(const AqueousMixture& mixture) -> unsigned
{
	return mixture.indicesNeutralSpecies().size();
}

auto localIndexChargedSpecies(const AqueousMixture& mixture, const std::string& name) -> Index
{
	const Index idx = indexSpecies(mixture, name);
	return find(idx, mixture.indicesChargedSpecies());
}

auto localIndexNeutralSpecies(const AqueousMixture& mixture, const std::string& name) -> Index
{
	const Index idx = indexSpecies(mixture, name);
	return find(idx, mixture.indicesNeutralSpecies());
}

auto localIndexCation(const AqueousMixture& mixture, const std::string& name) -> Index
{
	const Index idx = indexSpecies(mixture, name);
	return find(idx, indicesCations(mixture));
}

auto localIndexAnion(const AqueousMixture& mixture, const std::string& name) -> Index
{
	const Index idx = indexSpecies(mixture, name);
	return find(idx, indicesAnions(mixture));
}

auto indicesCations(const AqueousMixture& mixture) -> Indices
{
	Indices indices_cations;
	for(const auto& idx : mixture.indicesChargedSpecies())
		if(mixture.charges()[idx] > 0) indices_cations.push_back(idx);
	return indices_cations;
}

auto indicesAnions(const AqueousMixture& mixture) -> Indices
{
	Indices indices_anions;
	for(const auto& idx : mixture.indicesChargedSpecies())
		if(mixture.charges()[idx] < 0) indices_anions.push_back(idx);
	return indices_anions;
}

auto namesChargedSpecies(const AqueousMixture& mixture) -> std::vector<std::string>
{
	return extract(namesSpecies(mixture), mixture.indicesChargedSpecies());
}

auto namesNeutralSpecies(const AqueousMixture& mixture) -> std::vector<std::string>
{
	return extract(namesSpecies(mixture), mixture.indicesNeutralSpecies());
}

auto namesCations(const AqueousMixture& mixture) -> std::vector<std::string>
{
	return extract(namesSpecies(mixture), indicesCations(mixture));
}

auto namesAnions(const AqueousMixture& mixture) -> std::vector<std::string>
{
	return extract(namesSpecies(mixture), indicesAnions(mixture));
}

auto molalities(const AqueousMixture& mixture, const Vector& n) -> VectorResult
{
    const unsigned size = numSpecies(mixture);
    const Index index_water = mixture.indexWater();
    const double n_water = n[index_water];

	VectorResult res;
	res.func = 55.508 * n/n_water;
	res.grad = zeros(size, size);

    for(unsigned i = 0; i < size; ++i)
    	res.grad(i, i) = res.func[i]/n[i];

    for(unsigned i = 0; i < size; ++i)
    	res.grad(i, index_water) -= res.func[i]/n_water;

    return res;
}

auto molalitiesStoichiometric(const AqueousMixture& mixture, const VectorResult& m) -> VectorResult
{
	const Indices& indices_charged = mixture.indicesChargedSpecies();
	const Indices& indices_neutral = mixture.indicesNeutralSpecies();
	const Matrix& dissociation_mat = mixture.dissociationMatrix();
	
    // The molalities of the charged species
    VectorResult m_charged;
    m_charged.func = rows(indices_charged, m.func);
    m_charged.grad = rows(indices_charged, m.grad);

    // The molalities of the neutral species
    VectorResult m_neutral;
    m_neutral.func = rows(indices_neutral, m.func);
    m_neutral.grad = rows(indices_neutral, m.grad);

    // The stoichiometric molalities of the charged species
    VectorResult ms;
    ms.func = m_charged.func + dissociation_mat.t() * m_neutral.func;
    ms.grad = m_charged.grad + dissociation_mat.t() * m_neutral.grad;

    return ms;
}

auto ionicStrength(const AqueousMixture& mixture, const VectorResult& m) -> ScalarResult
{
    const unsigned num_species = numSpecies(mixture);
    const Vector& z = mixture.charges();
    ScalarResult res(num_species);
    res.func = 0.5 * arma::sum(z * z * m.func);
    for(unsigned j = 0; j < num_species; ++j)
        res.grad[j] = 0.5 * arma::sum(z * z * m.grad.col(j));
    return res;
}

auto ionicStrengthStoichiometric(const AqueousMixture& mixture, const VectorResult& ms) -> ScalarResult
{
    const unsigned num_species = numSpecies(mixture);
    const Vector& z_charged = mixture.chargesChargedSpecies();
    ScalarResult res(num_species);
    res.func = 0.5 * arma::sum(z_charged * z_charged * ms.func);
    for(unsigned j = 0; j < num_species; ++j)
        res.grad[j] = 0.5 * arma::sum(z_charged * z_charged * ms.grad.col(j));
    return res;
}

} // namespace Reaktor
