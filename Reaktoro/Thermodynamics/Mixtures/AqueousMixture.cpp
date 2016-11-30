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

#include "AqueousMixture.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Common/InterpolationUtils.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterThermoStateUtils.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterUtils.hpp>

namespace Reaktoro {
namespace internal {

auto indicesChargedSpecies(const std::vector<AqueousSpecies>& mixture) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < mixture.size(); ++i)
        if(mixture[i].charge() != 0)
            indices.push_back(i);
    return indices;
}

auto indicesNeutralSpecies(const std::vector<AqueousSpecies>& mixture) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < mixture.size(); ++i)
        if(mixture[i].charge() == 0)
            indices.push_back(i);
    return indices;
}

auto indicesCations(const std::vector<AqueousSpecies>& mixture) -> Indices
{
    Indices indices_cations;
    for(unsigned i = 0; i < mixture.size(); ++i)
        if(mixture[i].charge() > 0)
            indices_cations.push_back(i);
    return indices_cations;
}

auto indicesAnions(const std::vector<AqueousSpecies>& mixture) -> Indices
{
    Indices indices_anions;
    for(unsigned i = 0; i < mixture.size(); ++i)
        if(mixture[i].charge() < 0)
            indices_anions.push_back(i);
    return indices_anions;
}

auto dissociationMatrix(const std::vector<AqueousSpecies>& mixture) -> Matrix
{
    // The indices of the neutral and charged species
    const Indices indices_neutral = indicesNeutralSpecies(mixture);
    const Indices indices_charged = indicesChargedSpecies(mixture);

    // Gets the stoichiometry of the i-th charged species in the j-th neutral species
    auto stoichiometry = [&](unsigned i, unsigned j) -> double
    {
        const Index ineutral = indices_neutral[i];
        const Index icharged = indices_charged[j];
        const AqueousSpecies& neutral = mixture[ineutral];
        const AqueousSpecies& charged = mixture[icharged];
        const auto iter = neutral.dissociation().find(charged.name());
        return iter != neutral.dissociation().end() ? iter->second : 0.0;
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

auto defaultWaterDensityFunction() -> ThermoScalarFunction
{
    ThermoScalarFunction f = [=](double T, double P) -> ThermoScalar
    {
        return waterDensityWagnerPruss(T, P);
    };

    return f;
}

auto defaultWaterDielectricConstantFunction() -> ThermoScalarFunction
{
    WaterThermoState wts;
    WaterElectroState wes;

    ThermoScalarFunction f = [=](double T, double P) mutable -> ThermoScalar
    {
        wts = waterThermoStateHGK(T, P);
        wes = waterElectroStateJohnsonNorton(T, P, wts);
        return wes.epsilon;
    };

    return f;
}

} // namespace internal

AqueousMixture::AqueousMixture()
{}

AqueousMixture::AqueousMixture(const std::vector<AqueousSpecies>& species)
: GeneralMixture<AqueousSpecies>(species)
{
    // Initialize the index of the water species
    idx_water = indexSpeciesAny(alternativeWaterNames());

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

auto AqueousMixture::dissociationMatrix() const -> const Matrix&
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

auto AqueousMixture::chargesChargedSpecies() const -> Vector
{
    return rows(chargesSpecies(), indicesChargedSpecies());
}

auto AqueousMixture::chargesCations() const -> Vector
{
    return rows(chargesSpecies(), indicesCations());
}

auto AqueousMixture::chargesAnions() const -> Vector
{
    return rows(chargesSpecies(), indicesAnions());
}

auto AqueousMixture::molalities(const Vector& n) const -> ChemicalVector
{
    const unsigned num_species = numSpecies();

    // The molalities of the species and their partial derivatives
    ChemicalVector m(num_species);

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

auto AqueousMixture::stoichiometricMolalities(const ChemicalVector& m) const -> ChemicalVector
{
    // Auxiliary variables
    const unsigned num_species = numSpecies();
    const unsigned num_charged = numChargedSpecies();
    const unsigned num_neutral = numNeutralSpecies();

    // The molalities of the charged species
    ChemicalVector mc(num_charged, num_species);
    mc.val = rows(m.val, idx_charged_species);
    mc.ddn = rows(m.ddn, idx_charged_species);

    // The molalities of the neutral species
    ChemicalVector mn(num_neutral, num_species);
    mn.val = rows(m.val, idx_neutral_species);
    mn.ddn = rows(m.ddn, idx_neutral_species);

    // The stoichiometric molalities of the charged species
    ChemicalVector ms(num_charged, num_species);
    ms.val = mc.val + tr(dissociation_matrix) * mn.val;
    ms.ddn = mc.ddn + tr(dissociation_matrix) * mn.ddn;

    return ms;
}

auto AqueousMixture::effectiveIonicStrength(const ChemicalVector& m) const -> ChemicalScalar
{
    const unsigned num_species = numSpecies();
    const Vector z = chargesSpecies();

    ChemicalScalar Ie(num_species);
    Ie.val = 0.5 * sum(z % z % m.val);
    for(unsigned i = 0; i < num_species; ++i)
        Ie.ddn[i] = 0.5 * sum(z % z % m.ddn.col(i));

    return Ie;
}

auto AqueousMixture::stoichiometricIonicStrength(const ChemicalVector& ms) const -> ChemicalScalar
{
    const unsigned num_species = numSpecies();
    const Vector zc = chargesChargedSpecies();

    ChemicalScalar Is(num_species);
    Is.val = 0.5 * sum(zc % zc % ms.val);
    for(unsigned i = 0; i < num_species; ++i)
        Is.ddn[i] = 0.5 * sum(zc % zc % ms.ddn.col(i));

    return Is;
}

auto AqueousMixture::state(double T, double P, const Vector& n) const -> AqueousMixtureState
{
    AqueousMixtureState res;
    res.T = T;
    res.P = P;
    res.x = molarFractions(n);
    res.rho = rho(T, P);
    res.epsilon = epsilon(T, P);
    res.m  = molalities(n);
    res.ms = stoichiometricMolalities(res.m);
    res.Ie = effectiveIonicStrength(res.m);
    res.Is = stoichiometricIonicStrength(res.ms);
    return res;
}

} // namespace Reaktoro
