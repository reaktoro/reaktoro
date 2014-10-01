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

auto operator==(const SolutionState& l, const SolutionState& r) -> bool
{
    return l.T == r.T and l.P == r.P and arma::all(l.n == r.n);
}

auto waterIndex(const AqueousSolution& solution) -> Index
{
    auto comparer = [](const AqueousSpecies& species)
    {
        return species.name == "H2O(l)";
    };

    return std::find_if(solution.begin(), solution.end(), comparer) - solution.begin();
}

auto chargedSpeciesIndices(const AqueousSolution& solution) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < numSpecies(solution); ++i)
        if(solution[i].charge != 0)
            indices.push_back(i);
    return indices;
}

auto neutralSpeciesIndices(const AqueousSolution& solution) -> Indices
{
    Indices indices;
    for(unsigned i = 0; i < numSpecies(solution); ++i)
        if(solution[i].charge == 0)
            indices.push_back(i);
    return indices;
}

auto cationIndices(const AqueousSolution& solution) -> Indices
{
    Indices indices_cations;
    for(unsigned i = 0; i < numSpecies(solution); ++i)
        if(solution[i].charge > 0)
            indices_cations.push_back(i);
    return indices_cations;
}

auto anionIndices(const AqueousSolution& solution) -> Indices
{
    Indices indices_anions;
    for(unsigned i = 0; i < numSpecies(solution); ++i)
        if(solution[i].charge < 0)
            indices_anions.push_back(i);
    return indices_anions;
}

auto chargedSpeciesLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return find(idx, chargedSpeciesIndices(solution));
}

auto neutralSpeciesLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return find(idx, neutralSpeciesIndices(solution));
}

auto cationLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return find(idx, cationIndices(solution));
}

auto anionLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return find(idx, anionIndices(solution));
}

auto speciesNames(const AqueousSolution& solution) -> std::vector<std::string>
{
    const unsigned nspecies = numSpecies(solution);
    std::vector<std::string> names(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        names[i] = solution[i].name;
    return names;
}

auto chargedSpeciesNames(const AqueousSolution& solution) -> std::vector<std::string>
{
    return extract(speciesNames(solution), chargedSpeciesIndices(solution));
}

auto neutralSpeciesNames(const AqueousSolution& solution) -> std::vector<std::string>
{
    return extract(speciesNames(solution), neutralSpeciesIndices(solution));
}

auto cationNames(const AqueousSolution& solution) -> std::vector<std::string>
{
    return extract(speciesNames(solution), cationIndices(solution));
}

auto anionNames(const AqueousSolution& solution) -> std::vector<std::string>
{
    return extract(speciesNames(solution), anionIndices(solution));
}

auto chargedSpeciesCharges(const AqueousSolution& solution) -> Vector
{
    const arma::uvec icharged = chargedSpeciesIndices(solution);
    return speciesCharges(solution).rows(icharged);
}

auto cationCharges(const AqueousSolution& solution) -> Vector
{
    const arma::uvec ications = cationIndices(solution);
    return speciesCharges(solution).rows(ications);
}

auto anionCharges(const AqueousSolution& solution) -> Vector
{
    const arma::uvec ianions = anionIndices(solution);
    return speciesCharges(solution).rows(ianions);
}

auto dissociationMatrix(const AqueousSolution& solution) -> Matrix
{
    // The indices of the neutral and charged species
    const Indices indices_neutral = neutralSpeciesIndices(solution);
    const Indices indices_charged = chargedSpeciesIndices(solution);

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

auto updateSolutionState(SolutionState& state, double T, double P, const Vector& n) -> void
{
    state.T = T;
    state.P = P;
    state.n = n;
    state.x = molarFractions(n);
}

auto aqueousSolutionStateFunction(const AqueousSolution& solution) -> AqueousSolutionStateFunction
{
    // Auxiliary variables
    const Index    iH2O             = waterIndex(solution);
    const Vector   z                = speciesCharges(solution);
    const Vector   z_charged        = chargedSpeciesCharges(solution);
    const Indices  indices_charged  = chargedSpeciesIndices(solution);
    const Indices  indices_neutral  = neutralSpeciesIndices(solution);
    const unsigned num_species      = numSpecies(solution);
    const Matrix   dissociation_mat = dissociationMatrix(solution);

    // Auxiliary variables for performance reasons
    Vector m_charged, m_neutral;
    Matrix dmdn_charged, dmdn_neutral;

    const Vector zero = zeros(num_species);

    double Ie = 0, Is = 0;
    Vector dIedn, dIsdn;

    Vector m, ms;
    Matrix dmdn, dmsdn;

    // The state of the solution and some references for clarity reasons
    AqueousSolutionState state;

    AqueousSolutionStateFunction func = [=](double T, double P, const Vector& n) mutable -> AqueousSolutionState
    {
        // Update the temperature, pressure and amounts of the species
        updateSolutionState(state, T, P, n);

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

        // Computing the effective ionic strength of the aqueous solution
        Ie = 0.5 * arma::sum(z * z * m);
        for(unsigned i = 0; i < num_species; ++i)
            dIedn[i] = 0.5 * arma::sum(z * z * dmdn.col(i));

        // Computing the stoichiometric ionic strength of the aqueous solution
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

auto gaseousSolutionStateFunction(const GaseousSolution& solution) -> GaseousSolutionStateFunction
{
    GaseousSolutionState state;
    GaseousSolutionStateFunction func = [=](double T, double P, const Vector& n) mutable -> GaseousSolutionState
    {
        updateSolutionState(state, T, P, n);
        return state;
    };
    return func;
}

auto mineralSolutionStateFunction(const GaseousSolution& solution) -> MineralSolutionStateFunction
{
    MineralSolutionState st;
    MineralSolutionStateFunction func = [=](double T, double P, const Vector& n) mutable -> MineralSolutionState
    {
        updateSolutionState(st, T, P, n);
        return st;
    };
    return func;
}

} // namespace Reaktor
