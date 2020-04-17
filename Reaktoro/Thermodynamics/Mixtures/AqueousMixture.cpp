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

#include "AqueousMixture.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/InterpolationUtils.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Singletons/DissociationReactions.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>

namespace Reaktoro {
namespace internal {

auto defaultWaterDensityFunction() -> std::function<real(real, real)>
{
    const auto T = 298.15;
    const auto P = 1.0e5;
    const auto rho = waterLiquidDensityWagnerPruss(T, P);
    return [=](real T, real P) { return rho; };
}

auto defaultWaterDielectricConstantFunction() -> std::function<real(real, real)>
{
    const auto T = 298.15;
    const auto P = 1.0e5;
    const auto wts = waterThermoStateHGK(T, P, StateOfMatter::Liquid);
    const auto wes = waterElectroStateJohnsonNorton(T, P, wts);
    const auto epsilon = wes.epsilon;
    return [=](real T, real P) { return epsilon; };
}

} // namespace internal

AqueousMixture::AqueousMixture()
{}

AqueousMixture::AqueousMixture(const std::vector<Species>& species)
: GeneralMixture(species)
{
    // Initialize the index related data
    initializeIndices(species);

    // Initialize the dissociation matrix of the neutral species w.r.t. the charged species
    initializeDissociationMatrix(species);

    // Initialize the density function for water
    rho = rho_default = internal::defaultWaterDensityFunction();

    // Initialize the dielectric constant function for water
    epsilon = epsilon_default = internal::defaultWaterDielectricConstantFunction();
}

AqueousMixture::~AqueousMixture()
{}

auto AqueousMixture::setWaterDensity(std::function<real(real, real)> rho) -> void
{
    this->rho = rho;
}

auto AqueousMixture::setWaterDielectricConstant(std::function<real(real, real)> epsilon) -> void
{
    this->epsilon = epsilon;
}

auto AqueousMixture::setInterpolationPoints(const std::vector<real>& temperatures, const std::vector<real>& pressures) -> void
{
    rho = interpolate(temperatures, pressures, rho_default);
    epsilon = interpolate(temperatures, pressures, epsilon_default);
}

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

auto AqueousMixture::dissociationMatrix() const -> MatrixXdConstRef
{
    return dissociation_matrix;
}

auto AqueousMixture::indexNeutralSpecies(std::string name) const -> Index
{
    const Index idx = indexSpecies(name);
    return index(idx_neutral_species, idx);
}

auto AqueousMixture::indexNeutralSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    const Index idx = indexSpeciesAny(names);
    return index(idx_neutral_species, idx);
}

auto AqueousMixture::indexChargedSpecies(std::string name) const -> Index
{
    const Index idx = indexSpecies(name);
    return index(idx_charged_species, idx);
}

auto AqueousMixture::indexChargedSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    const Index idx = indexSpeciesAny(names);
    return index(idx_charged_species, idx);
}

auto AqueousMixture::indexCation(std::string name) const -> Index
{
    const Index idx = indexSpecies(name);
    return index(idx_cations, idx);
}

auto AqueousMixture::indexAnion(std::string name) const -> Index
{
    const Index idx = indexSpecies(name);
    return index(idx_anions, idx);
}

auto AqueousMixture::namesNeutralSpecies() const -> std::vector<std::string>
{
    return extract(namesSpecies(), indicesNeutralSpecies());
}

auto AqueousMixture::namesChargedSpecies() const -> std::vector<std::string>
{
    return extract(namesSpecies(), indicesChargedSpecies());
}

auto AqueousMixture::namesCations() const -> std::vector<std::string>
{
    return extract(namesSpecies(), indicesCations());
}

auto AqueousMixture::namesAnions() const -> std::vector<std::string>
{
    return extract(namesSpecies(), indicesAnions());
}

auto AqueousMixture::chargesChargedSpecies() const -> ArrayXr
{
    return charges()(indicesChargedSpecies());
}

auto AqueousMixture::chargesCations() const -> ArrayXr
{
    return charges()(indicesCations());
}

auto AqueousMixture::chargesAnions() const -> ArrayXr
{
    return charges()(indicesAnions());
}

auto AqueousMixture::molalities(ArrayXrConstRef x) const -> ArrayXr
{
    const auto xw = x[idx_water];
    if(xw == 0.0)
        return ArrayXr::Zero(x.size());
    return x/(waterMolarMass * xw);
}

auto AqueousMixture::stoichiometricMolalities(ArrayXrConstRef m) const -> ArrayXr
{
    // The molalities of the charged species
    const auto mc = m(idx_charged_species).matrix(); // convert from array to matrix expression

    // The molalities of the neutral species
    const auto mn = m(idx_neutral_species).matrix(); // convert from array to matrix expression

    // The stoichiometric molalities of the charged species
    const auto ms = mc + dissociation_matrix.transpose() * mn;

    return ms;
}

auto AqueousMixture::effectiveIonicStrength(ArrayXrConstRef m) const -> real
{
    const auto num_species = numSpecies();
    const auto z = charges();
    return 0.5 * (z * z * m).sum();
}

auto AqueousMixture::stoichiometricIonicStrength(ArrayXrConstRef ms) const -> real
{
    const auto num_species = numSpecies();
    const auto zc = chargesChargedSpecies();
    return 0.5 * (zc * zc * ms).sum();
}

auto AqueousMixture::state(real T, real P, ArrayXrConstRef x) const -> AqueousMixtureState
{
    AqueousMixtureState res;
    res.T = T;
    res.P = P;
    res.rho = rho(T, P);
    res.epsilon = epsilon(T, P);
    res.m  = molalities(x);
    res.ms = stoichiometricMolalities(res.m);
    res.Ie = effectiveIonicStrength(res.m);
    res.Is = stoichiometricIonicStrength(res.ms);
    return res;
}

auto AqueousMixture::initializeIndices(const std::vector<Species>& species) -> void
{
    // Initialize the index of the water species
    idx_water = indexSpeciesAny(alternativeWaterNames());

    // Ensure water species is present
    Assert(idx_water < numSpecies(),
        "Could not initialize the aqueous mixture.",
        "You probably forgot to add water species in the definition of the aqueous phase (e.g., H2O(l)).");

    // Initialize the indices of the charged and neutral species
    for(unsigned i = 0; i < species.size(); ++i)
    {
        if(i == idx_water) continue; // Skip if water species
        if(species[i].charge() == 0)
            idx_neutral_species.push_back(i); // Current species is neutral
        else
        {
            idx_charged_species.push_back(i); // Current species is charged
            if(species[i].charge() > 0) idx_cations.push_back(i); // Current species is a cation
            else idx_anions.push_back(i); // Current species is an anion
        }
    }
}

auto AqueousMixture::initializeDissociationMatrix(const std::vector<Species>& species) -> void
{
    // Return the stoichiometry of the i-th charged species in the j-th neutral species
    auto stoichiometry = [&](Index i, Index j) -> real
    {
        const auto ineutral = idx_neutral_species[i];
        const auto icharged = idx_charged_species[j];
        const auto neutral = species[ineutral].formula();
        const auto charged = species[icharged].formula();
        return DissociationReactions::coefficient(neutral, charged); // TODO: Performance can be improved here by avoiding recreation of ChemicalFormula objects.
    };

    // Assemble the dissociation matrix of the neutral species with respect to the charged species
    const Index num_charged_species = idx_charged_species.size();
    const Index num_neutral_species = idx_neutral_species.size();
    dissociation_matrix.resize(num_neutral_species, num_charged_species);
    for(Index i = 0; i < num_neutral_species; ++i)
        for(Index j = 0; j < num_charged_species; ++j)
            dissociation_matrix(i, j) = stoichiometry(i, j);
}

} // namespace Reaktoro
