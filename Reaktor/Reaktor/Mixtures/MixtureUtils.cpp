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

#include "MixtureUtils.hpp"

// Reaktor includes
#include <Reaktor/Common/MatrixUtils.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>

namespace Reaktor {

auto indexWater(const AqueousMixture& mixture) -> Index
{
    auto comparer = [](const AqueousSpecies& species)
    {
        return species.name == "H2O(l)";
    };

    return std::find_if(mixture.begin(), mixture.end(), comparer) - mixture.begin();
}

auto indicesChargedSpecies(const AqueousMixture& mixture) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].charge != 0)
            indices.push_back(i);
    return indices;
}

auto indicesNeutralSpecies(const AqueousMixture& mixture) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].charge == 0)
            indices.push_back(i);
    return indices;
}

auto indicesCations(const AqueousMixture& mixture) -> Indices
{
    Indices indices_cations;
    for(unsigned i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].charge > 0)
            indices_cations.push_back(i);
    return indices_cations;
}

auto indicesAnions(const AqueousMixture& mixture) -> Indices
{
    Indices indices_anions;
    for(unsigned i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].charge < 0)
            indices_anions.push_back(i);
    return indices_anions;
}

auto localIndexChargedSpecies(const AqueousMixture& mixture, const std::string& name) -> Index
{
    const Index idx = indexSpecies(mixture, name);
    return find(idx, indicesChargedSpecies(mixture));
}

auto localIndexNeutralSpecies(const AqueousMixture& mixture, const std::string& name) -> Index
{
    const Index idx = indexSpecies(mixture, name);
    return find(idx, indicesNeutralSpecies(mixture));
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

auto namesSpecies(const AqueousMixture& mixture) -> std::vector<std::string>
{
    const unsigned nspecies = numSpecies(mixture);
    std::vector<std::string> names(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        names[i] = mixture[i].name;
    return names;
}

auto namesChargedSpecies(const AqueousMixture& mixture) -> std::vector<std::string>
{
    return extract(namesSpecies(mixture), indicesChargedSpecies(mixture));
}

auto namesNeutralSpecies(const AqueousMixture& mixture) -> std::vector<std::string>
{
    return extract(namesSpecies(mixture), indicesNeutralSpecies(mixture));
}

auto namesCations(const AqueousMixture& mixture) -> std::vector<std::string>
{
    return extract(namesSpecies(mixture), indicesCations(mixture));
}

auto namesAnions(const AqueousMixture& mixture) -> std::vector<std::string>
{
    return extract(namesSpecies(mixture), indicesAnions(mixture));
}

auto chargesChargedSpecies(const AqueousMixture& mixture) -> Vector
{
    const arma::uvec icharged = indicesChargedSpecies(mixture);
    return chargesSpecies(mixture).rows(icharged);
}

auto chargesCations(const AqueousMixture& mixture) -> Vector
{
    const arma::uvec ications = indicesCations(mixture);
    return chargesSpecies(mixture).rows(ications);
}

auto chargesAnions(const AqueousMixture& mixture) -> Vector
{
    const arma::uvec ianions = indicesAnions(mixture);
    return chargesSpecies(mixture).rows(ianions);
}

auto dissociationMatrix(const AqueousMixture& mixture) -> Matrix
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
        const auto iter = neutral.dissociation.find(charged.name);
        return iter != neutral.dissociation.end() ? iter->second : 0.0;
    };

    // Assemble the dissociation matrix of the neutral species with respect to the charged species
    const unsigned num_charged_species = indices_charged.size();
    const unsigned num_neutral_species = indices_neutral.size();
    Matrix dissociation_matrix = zeros(num_charged_species, num_neutral_species);
    for (unsigned i = 0; i < num_neutral_species; ++i)
        for (unsigned j = 0; j < num_charged_species; ++j)
            dissociation_matrix(i, j) = stoichiometry(i, j);
    return dissociation_matrix;
}

auto molarFractions(const Vector& n) -> ThermoVector
{
    const unsigned nspecies = n.size();
    const double nt = arma::sum(n);
    ThermoVector x = ThermoVector::zero(nspecies, nspecies);
    x.val = n/nt;
    for(unsigned i = 0; i < nspecies; ++i)
    {
        x.ddn.row(i).fill(-x.val[i]/nt);
        x.ddn(i, i) += 1.0/nt;
    }
    return x;
}

auto molalities(const AqueousMixture& mixture, const Vector& n) -> ThermoVector
{
    const unsigned size = numSpecies(mixture);
    const Index index_water = indexWater(mixture);
    const double n_water = n[index_water];

    ThermoVector res;
    res.val = 55.508 * n/n_water;
    res.ddn = zeros(size, size);

    for(unsigned i = 0; i < size; ++i)
        res.ddn(i, i) = res.val[i]/n[i];

    for(unsigned i = 0; i < size; ++i)
        res.ddn(i, index_water) -= res.val[i]/n_water;

    return res;
}

auto molalitiesStoichiometric(const AqueousMixture& mixture, const ThermoVector& m) -> ThermoVector
{
    const Indices& indices_charged = indicesChargedSpecies(mixture);
    const Indices& indices_neutral = indicesNeutralSpecies(mixture);
    const Matrix& dissociation_mat = dissociationMatrix(mixture);

    // The molalities of the charged species
    ThermoVector m_charged;
    m_charged.val = rows(indices_charged, m.val);
    m_charged.ddn = rows(indices_charged, m.ddn);

    // The molalities of the neutral species
    ThermoVector m_neutral;
    m_neutral.val = rows(indices_neutral, m.val);
    m_neutral.ddn = rows(indices_neutral, m.ddn);

    // The stoichiometric molalities of the charged species
    ThermoVector ms;
    ms.val = m_charged.val + dissociation_mat.t() * m_neutral.val;
    ms.ddn = m_charged.ddn + dissociation_mat.t() * m_neutral.ddn;

    return ms;
}

auto ionicStrength(const AqueousMixture& mixture, const ThermoVector& m) -> ThermoScalar
{
    const unsigned num_species = numSpecies(mixture);
    const Vector& z = chargesSpecies(mixture);
    ThermoScalar res(num_species);
    res.val = 0.5 * arma::sum(z * z * m.val);
    for(unsigned j = 0; j < num_species; ++j)
        res.ddn[j] = 0.5 * arma::sum(z * z * m.ddn.col(j));
    return res;
}

auto ionicStrengthStoichiometric(const AqueousMixture& mixture, const ThermoVector& ms) -> ThermoScalar
{
    const unsigned num_species = numSpecies(mixture);
    const Vector& z_charged = chargesChargedSpecies(mixture);
    ThermoScalar res(num_species);
    res.val = 0.5 * arma::sum(z_charged * z_charged * ms.val);
    for(unsigned j = 0; j < num_species; ++j)
        res.ddn[j] = 0.5 * arma::sum(z_charged * z_charged * ms.ddn.col(j));
    return res;
}

auto operator==(const MixtureState& l, const MixtureState& r) -> bool
{
    return l.T == r.T and l.P == r.P and arma::all(l.n == r.n);
}

auto updateMixtureState(MixtureState& state, double T, double P, const Vector& n) -> void
{
    state.T = T;
    state.P = P;
    state.n = n;
    state.x = molarFractions(n);
}

auto aqueousMixtureStateFunction(const AqueousMixture& mixture) -> AqueousMixtureStateFunction
{
    // Auxiliary variables
    const Index    iH2O             = indexWater(mixture);
    const Vector   z                = chargesSpecies(mixture);
    const Vector   z_charged        = chargesChargedSpecies(mixture);
    const Indices  indices_charged  = indicesChargedSpecies(mixture);
    const Indices  indices_neutral  = indicesNeutralSpecies(mixture);
    const unsigned num_species      = numSpecies(mixture);
    const unsigned num_charged      = indices_charged.size();
    const Matrix   dissociation_mat = dissociationMatrix(mixture);

    // Auxiliary variables for performance reasons
    ThermoVector m_charged;
    ThermoVector m_neutral;

    // The state of the mixture and some references for clarity reasons
    AqueousMixtureState state;

    AqueousMixtureStateFunction func = [=](double T, double P, const Vector& n) mutable -> AqueousMixtureState
    {
        // Update the temperature, pressure and amounts of the species
        updateMixtureState(state, T, P, n);

        // Computing the molalities of the species
        state.m = ThermoVector::zero(num_species, num_species);
        state.m.val = 55.508 * n/n[iH2O];
        for(unsigned i = 0; i < num_species; ++i)
        {
            state.m.ddn(i, i)     = state.m.val[i]/n[i];
            state.m.ddn(i, iH2O) -= state.m.val[i]/n[iH2O];
        }

        // Computing the stoichiometric molalities of the charged species
        state.ms = ThermoVector::zero(num_charged, num_species);
        m_charged.val = rows(indices_charged, state.m.val);
        m_charged.ddn = rows(indices_charged, state.m.ddn);
        m_neutral.val = rows(indices_neutral, state.m.val);
        m_neutral.ddn = rows(indices_neutral, state.m.ddn);
        state.ms.val = m_charged.val + dissociation_mat.t() * m_neutral.val;
        state.ms.ddn = m_charged.ddn + dissociation_mat.t() * m_neutral.ddn;

        // Computing the effective ionic strength of the aqueous mixture
        state.Ie = ThermoScalar::zero(num_species);
        state.Ie.val = 0.5 * arma::sum(z * z * state.m.val);
        for(unsigned i = 0; i < num_species; ++i)
            state.Ie.ddn[i] = 0.5 * arma::sum(z * z * state.m.ddn.col(i));

        // Computing the stoichiometric ionic strength of the aqueous mixture
        state.Is = ThermoScalar::zero(num_species);
        state.Is.val = 0.5 * arma::sum(z_charged * z_charged * state.ms.val);
        for(unsigned i = 0; i < num_species; ++i)
            state.Is.ddn[i] = 0.5 * arma::sum(z_charged * z_charged * state.ms.ddn.col(i));

        return state;
    };

    return func;
}

auto gaseousMixtureStateFunction(const GaseousMixture& mixture) -> GaseousMixtureStateFunction
{
    GaseousMixtureState state;
    GaseousMixtureStateFunction func = [=](double T, double P, const Vector& n) mutable -> GaseousMixtureState
    {
        updateMixtureState(state, T, P, n);
        return state;
    };
    return func;
}

auto mineralMixtureStateFunction(const GaseousMixture& mixture) -> MineralMixtureStateFunction
{
    MineralMixtureState st;
    MineralMixtureStateFunction func = [=](double T, double P, const Vector& n) mutable -> MineralMixtureState
    {
        updateMixtureState(st, T, P, n);
        return st;
    };
    return func;
}

} // namespace Reaktor
