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

#include "ActivityUtils.hpp"

// Reaktor includes
#include <Reaktor/Common/MatrixUtils.hpp>
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>

namespace Reaktor {

auto operator==(const MixtureState& l, const MixtureState& r) -> bool
{
    return l.T == r.T and l.P == r.P and arma::all(l.n == r.n);
}

auto waterIndex(const AqueousMixture& mixture) -> Index
{
    auto comparer = [](const AqueousSpecies& species)
    {
        return species.name == "H2O(l)";
    };

    return std::find_if(mixture.begin(), mixture.end(), comparer) - mixture.begin();
}

auto chargedSpeciesIndices(const AqueousMixture& mixture) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].charge != 0)
            indices.push_back(i);
    return indices;
}

auto neutralSpeciesIndices(const AqueousMixture& mixture) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].charge == 0)
            indices.push_back(i);
    return indices;
}

auto cationIndices(const AqueousMixture& mixture) -> Indices
{
    Indices indices_cations;
    for(unsigned i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].charge > 0)
            indices_cations.push_back(i);
    return indices_cations;
}

auto anionIndices(const AqueousMixture& mixture) -> Indices
{
    Indices indices_anions;
    for(unsigned i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].charge < 0)
            indices_anions.push_back(i);
    return indices_anions;
}

auto chargedSpeciesLocalIndex(const AqueousMixture& mixture, const std::string& name) -> Index
{
    const Index idx = speciesIndex(mixture, name);
    return find(idx, chargedSpeciesIndices(mixture));
}

auto neutralSpeciesLocalIndex(const AqueousMixture& mixture, const std::string& name) -> Index
{
    const Index idx = speciesIndex(mixture, name);
    return find(idx, neutralSpeciesIndices(mixture));
}

auto cationLocalIndex(const AqueousMixture& mixture, const std::string& name) -> Index
{
    const Index idx = speciesIndex(mixture, name);
    return find(idx, cationIndices(mixture));
}

auto anionLocalIndex(const AqueousMixture& mixture, const std::string& name) -> Index
{
    const Index idx = speciesIndex(mixture, name);
    return find(idx, anionIndices(mixture));
}

auto speciesNames(const AqueousMixture& mixture) -> std::vector<std::string>
{
    const unsigned nspecies = numSpecies(mixture);
    std::vector<std::string> names(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        names[i] = mixture[i].name;
    return names;
}

auto chargedSpeciesNames(const AqueousMixture& mixture) -> std::vector<std::string>
{
    return extract(speciesNames(mixture), chargedSpeciesIndices(mixture));
}

auto neutralSpeciesNames(const AqueousMixture& mixture) -> std::vector<std::string>
{
    return extract(speciesNames(mixture), neutralSpeciesIndices(mixture));
}

auto cationNames(const AqueousMixture& mixture) -> std::vector<std::string>
{
    return extract(speciesNames(mixture), cationIndices(mixture));
}

auto anionNames(const AqueousMixture& mixture) -> std::vector<std::string>
{
    return extract(speciesNames(mixture), anionIndices(mixture));
}

auto chargedSpeciesCharges(const AqueousMixture& mixture) -> Vector
{
    const arma::uvec icharged = chargedSpeciesIndices(mixture);
    return speciesCharges(mixture).rows(icharged);
}

auto cationCharges(const AqueousMixture& mixture) -> Vector
{
    const arma::uvec ications = cationIndices(mixture);
    return speciesCharges(mixture).rows(ications);
}

auto anionCharges(const AqueousMixture& mixture) -> Vector
{
    const arma::uvec ianions = anionIndices(mixture);
    return speciesCharges(mixture).rows(ianions);
}

auto dissociationMatrix(const AqueousMixture& mixture) -> Matrix
{
    // The indices of the neutral and charged species
    const Indices indices_neutral = neutralSpeciesIndices(mixture);
    const Indices indices_charged = chargedSpeciesIndices(mixture);

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

auto molarFractions(const Vector& n) -> ChemicalVector
{
    const unsigned nspecies = n.size();
    const double nt = arma::sum(n);
    Vector x = n/nt;
    Matrix dxdt = zeros(nspecies);
    Matrix dxdp = zeros(nspecies);
    Matrix dxdn = zeros(nspecies, nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
    {
        dxdn.row(i).fill(-x[i]/nt);
        dxdn(i, i) += 1.0/nt;
    }
    return {x, dxdt, dxdp, dxdn};
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
    const Index    iH2O             = waterIndex(mixture);
    const Vector   z                = speciesCharges(mixture);
    const Vector   z_charged        = chargedSpeciesCharges(mixture);
    const Indices  indices_charged  = chargedSpeciesIndices(mixture);
    const Indices  indices_neutral  = neutralSpeciesIndices(mixture);
    const unsigned num_species      = numSpecies(mixture);
    const Matrix   dissociation_mat = dissociationMatrix(mixture);

    // Auxiliary variables for performance reasons
    Vector m_charged, m_neutral;
    Matrix dmdn_charged, dmdn_neutral;

    const Vector zero = zeros(num_species);

    double Ie = 0, Is = 0;
    Vector dIedn, dIsdn;

    Vector m, ms;
    Matrix dmdn, dmsdn;

    // The state of the mixture and some references for clarity reasons
    AqueousMixtureState state;

    AqueousMixtureStateFunction func = [=](double T, double P, const Vector& n) mutable -> AqueousMixtureState
    {
        // Update the temperature, pressure and amounts of the species
        updateMixtureState(state, T, P, n);

        // Computing the molalities of the species
        m = 55.508 * n/n[iH2O];
        dmdn = zeros(num_species, num_species);
        for(unsigned i = 0; i < num_species; ++i)
        {
            dmdn(i, i)     = m[i]/n[i];
            dmdn(i, iH2O) -= m[i]/n[iH2O];
        }

        // Computing the stoichiometric molalities of the charged species
        m_charged = rows(indices_charged, m);
        m_neutral = rows(indices_neutral, m);
        dmdn_charged = rows(indices_charged, dmdn);
        dmdn_neutral = rows(indices_neutral, dmdn);

        ms = m_charged + dissociation_mat.t() * m_neutral;
        dmsdn = dmdn_charged + dissociation_mat.t() * dmdn_neutral;

        // Computing the effective ionic strength of the aqueous mixture
        Ie = 0.5 * arma::sum(z * z * m);
        for(unsigned i = 0; i < num_species; ++i)
            dIedn[i] = 0.5 * arma::sum(z * z * dmdn.col(i));

        // Computing the stoichiometric ionic strength of the aqueous mixture
        Is = 0.5 * arma::sum(z_charged * z_charged * ms);
        for(unsigned i = 0; i < num_species; ++i)
            dIsdn[i] = 0.5 * arma::sum(z_charged * z_charged * dmsdn.col(i));

        state.m  = ChemicalVector(m, zero, zero, dmdn);
        state.ms = ChemicalVector(ms, zero, zero, dmsdn);
        state.Ie = ChemicalScalar(Ie, 0.0, 0.0, dIedn);
        state.Is = ChemicalScalar(Is, 0.0, 0.0, dIsdn);

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
