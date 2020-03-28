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

#include "AqueousMixture.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/InterpolationUtils.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>

namespace Reaktoro {
namespace internal {

auto defaultWaterDensityFunction() -> ThermoScalarFunction
{
    const auto T = 298.15;
    const auto P = 1.0e5;
    const auto rho = waterLiquidDensityWagnerPruss(T, P);
    return [=](double T, double P) { return rho; };
}

auto defaultWaterDielectricConstantFunction() -> ThermoScalarFunction
{
    const auto T = 298.15;
    const auto P = 1.0e5;
    const auto wts = waterThermoStateHGK(T, P, StateOfMatter::Liquid);
    const auto wes = waterElectroStateJohnsonNorton(T, P, wts);
    const auto epsilon = wes.epsilon;
    return [=](double T, double P) { return epsilon; };
}

} // namespace internal

AqueousMixture::AqueousMixture()
{}

AqueousMixture::AqueousMixture(const std::vector<AqueousSpecies>& species)
: GeneralMixture<AqueousSpecies>(species)
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

auto AqueousMixture::setWaterDensity(const ThermoScalarFunction& rho) -> void
{
    this->rho = rho;
}

auto AqueousMixture::setWaterDielectricConstant(const ThermoScalarFunction& epsilon) -> void
{
    this->epsilon = epsilon;
}

auto AqueousMixture::setInterpolationPoints(const std::vector<double>& temperatures, const std::vector<double>& pressures) -> void
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
    return index(idx, idx_neutral_species);
}

auto AqueousMixture::indexNeutralSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    const Index idx = indexSpeciesAny(names);
    return index(idx, idx_neutral_species);
}

auto AqueousMixture::indexChargedSpecies(std::string name) const -> Index
{
    const Index idx = indexSpecies(name);
    return index(idx, idx_charged_species);
}

auto AqueousMixture::indexChargedSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    const Index idx = indexSpeciesAny(names);
    return index(idx, idx_charged_species);
}

auto AqueousMixture::indexCation(std::string name) const -> Index
{
    const Index idx = indexSpecies(name);
    return index(idx, idx_cations);
}

auto AqueousMixture::indexAnion(std::string name) const -> Index
{
    const Index idx = indexSpecies(name);
    return index(idx, idx_anions);
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

auto AqueousMixture::chargesChargedSpecies() const -> VectorXr
{
    return rows(chargesSpecies(), indicesChargedSpecies());
}

auto AqueousMixture::chargesCations() const -> VectorXr
{
    return rows(chargesSpecies(), indicesCations());
}

auto AqueousMixture::chargesAnions() const -> VectorXr
{
    return rows(chargesSpecies(), indicesAnions());
}

auto AqueousMixture::molalities(VectorXrConstRef n) const -> VectorXd
{
    const unsigned num_species = numSpecies();

    // The molalities of the species and their partial derivatives
    VectorXd m(num_species);

    // The molar amount of water
    const double nw = n[idx_water];

    // Check if the molar amount of water is zero
    if(nw == 0.0)
        return m;

    const double kgH2O = nw * waterMolarMass;

    m.val = n/kgH2O;
    for(unsigned i = 0; i < num_species; ++i)
    {
        m.ddn(i, i) = 1.0/kgH2O;
        m.ddn(i, idx_water) -= m.val[i]/nw;
    }

    return m;
}

auto AqueousMixture::stoichiometricMolalities(const VectorXd& m) const -> VectorXd
{
    // Auxiliary variables
    const unsigned num_species = numSpecies();
    const unsigned num_charged = numChargedSpecies();
    const unsigned num_neutral = numNeutralSpecies();

    // The molalities of the charged species
    VectorXd mc(num_charged, num_species);
    mc.val = rows(m.val, idx_charged_species);
    mc.ddn = rows(m.ddn, idx_charged_species);

    // The molalities of the neutral species
    VectorXd mn(num_neutral, num_species);
    mn.val = rows(m.val, idx_neutral_species);
    mn.ddn = rows(m.ddn, idx_neutral_species);

    // The stoichiometric molalities of the charged species
    VectorXd ms(num_charged, num_species);
    ms.val = mc.val + tr(dissociation_matrix) * mn.val;
    ms.ddn = mc.ddn + tr(dissociation_matrix) * mn.ddn;

    return ms;
}

auto AqueousMixture::effectiveIonicStrength(const VectorXd& m) const -> real
{
    const unsigned num_species = numSpecies();
    const VectorXr z = chargesSpecies();

    real Ie(num_species);
    Ie.val = 0.5 * sum(z % z % m.val);
    for(unsigned i = 0; i < num_species; ++i)
        Ie.ddn[i] = 0.5 * sum(z % z % m.ddn.col(i));

    return Ie;
}

auto AqueousMixture::stoichiometricIonicStrength(const VectorXd& ms) const -> real
{
    const unsigned num_species = numSpecies();
    const VectorXr zc = chargesChargedSpecies();

    real Is(num_species);
    Is.val = 0.5 * sum(zc % zc % ms.val);
    for(unsigned i = 0; i < num_species; ++i)
        Is.ddn[i] = 0.5 * sum(zc % zc % ms.ddn.col(i));

    return Is;
}

auto AqueousMixture::state(Temperature T, Pressure P, VectorXrConstRef n) const -> AqueousMixtureState
{
    AqueousMixtureState res;
    res.T = T;
    res.P = P;
    res.x = moleFractions(n);
    res.rho = rho(T, P);
    res.epsilon = epsilon(T, P);
    res.m  = molalities(n);
    res.ms = stoichiometricMolalities(res.m);
    res.Ie = effectiveIonicStrength(res.m);
    res.Is = stoichiometricIonicStrength(res.ms);
    return res;
}

auto AqueousMixture::initializeIndices(const std::vector<AqueousSpecies>& species) -> void
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

auto AqueousMixture::initializeDissociationMatrix(const std::vector<AqueousSpecies>& species) -> void
{
    // Return the stoichiometry of the i-th charged species in the j-th neutral species
    auto stoichiometry = [&](Index i, Index j) -> double
    {
        const Index ineutral = idx_neutral_species[i];
        const Index icharged = idx_charged_species[j];
        const AqueousSpecies& neutral = species[ineutral];
        const AqueousSpecies& charged = species[icharged];
        const auto iter = neutral.dissociation().find(charged.name());
        return iter != neutral.dissociation().end() ? iter->second : 0.0;
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
