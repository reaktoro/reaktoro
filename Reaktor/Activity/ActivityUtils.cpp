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
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>
#include <Reaktor/Thermodynamics/WaterConstants.hpp>

namespace Reaktor {

auto operator==(const SolutionState& l, const SolutionState& r) -> bool
{
    return l.T == r.T and l.P == r.P and l.n == r.n;
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
    return index(idx, chargedSpeciesIndices(solution));
}

auto neutralSpeciesLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return index(idx, neutralSpeciesIndices(solution));
}

auto cationLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return index(idx, cationIndices(solution));
}

auto anionLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index
{
    const Index idx = speciesIndex(solution, name);
    return index(idx, anionIndices(solution));
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
    return rows(speciesCharges(solution), chargedSpeciesIndices(solution));
}

auto cationCharges(const AqueousSolution& solution) -> Vector
{
    return rows(speciesCharges(solution), cationIndices(solution));
}

auto anionCharges(const AqueousSolution& solution) -> Vector
{
    return rows(speciesCharges(solution), anionIndices(solution));
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
    Matrix dissociation_matrix = zeros(num_neutral_species, num_charged_species);
    for(unsigned i = 0; i < num_neutral_species; ++i)
        for(unsigned j = 0; j < num_charged_species; ++j)
            dissociation_matrix(i, j) = stoichiometry(i, j);
    return dissociation_matrix;
}

auto molarFractions(const Vector& n) -> ChemicalVector
{
    const unsigned nspecies = n.size();
    const double nt = n.sum();
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
    const auto iH2O                = waterIndex(solution);
    const auto z                   = speciesCharges(solution);
    const auto z_charged           = chargedSpeciesCharges(solution);
    const auto indices_charged     = chargedSpeciesIndices(solution);
    const auto indices_neutral     = neutralSpeciesIndices(solution);
    const auto num_species         = numSpecies(solution);
    const auto num_charged_species = z_charged.size();
    const auto dissociation_mat    = dissociationMatrix(solution);

    // Auxiliary variables for performance reasons
    Vector m_val_charged;
    Vector m_val_neutral;
    Matrix m_ddn_charged;
    Matrix m_ddn_neutral;

    double Ie_val = 0.0;
    double Ie_ddt = 0.0;
    double Ie_ddp = 0.0;
    Vector Ie_ddn = zeros(num_species);

    double Is_val = 0.0;
    double Is_ddt = 0.0;
    double Is_ddp = 0.0;
    Vector Is_ddn = zeros(num_species);

    Vector m_val = zeros(num_species);
    Vector m_ddt = zeros(num_species);
    Vector m_ddp = zeros(num_species);
    Matrix m_ddn = zeros(num_species, num_species);

    Vector ms_val = zeros(num_charged_species);
    Vector ms_ddt = zeros(num_charged_species);
    Vector ms_ddp = zeros(num_charged_species);
    Matrix ms_ddn = zeros(num_charged_species, num_species);

    // The state of the solution and some references for clarity reasons
    AqueousSolutionState state;

    AqueousSolutionStateFunction func = [=](double T, double P, const Vector& n) mutable -> AqueousSolutionState
    {
        // Update the temperature, pressure and amounts of the species
        updateSolutionState(state, T, P, n);

        // Computing the molalities of the species
        m_val = n/(n[iH2O] * waterMolarMass);
        m_ddn = zeros(num_species, num_species);
        for(unsigned i = 0; i < num_species; ++i)
        {
            m_ddn(i, i)     = m_val[i]/n[i];
            m_ddn(i, iH2O) -= m_val[i]/n[iH2O];
        }

        // Computing the stoichiometric molalities of the charged species
        m_val_charged = rows(m_val, indices_charged);
        m_val_neutral = rows(m_val, indices_neutral);
        m_ddn_charged = rows(m_ddn, indices_charged);
        m_ddn_neutral = rows(m_ddn, indices_neutral);

        ms_val = m_val_charged + dissociation_mat.transpose() * m_val_neutral;
        ms_ddn = m_ddn_charged + dissociation_mat.transpose() * m_ddn_neutral;

        // Computing the effective ionic strength of the aqueous solution
        Ie_val = 0.5 * (z.array() * z.array() * m_val.array()).sum();
        for(unsigned i = 0; i < num_species; ++i)
            Ie_ddn[i] = 0.5 * (z.array() * z.array() * m_ddn.col(i).array()).sum();

        // Computing the stoichiometric ionic strength of the aqueous solution
        Is_val = 0.5 * (z_charged.array() * z_charged.array() * ms_val.array()).sum();
        for(unsigned i = 0; i < num_species; ++i)
            Is_ddn[i] = 0.5 * (z_charged.array() * z_charged.array() * ms_ddn.col(i).array()).sum();

        state.Ie = ChemicalScalar(Ie_val, Ie_ddt, Ie_ddp, Ie_ddn);
        state.Is = ChemicalScalar(Is_val, Is_ddt, Is_ddp, Is_ddn);
        state.m  = ChemicalVector(m_val, m_ddt, m_ddp, m_ddn);
        state.ms = ChemicalVector(ms_val, ms_ddt, ms_ddp, ms_ddn);

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

auto mineralSolutionStateFunction(const MineralSolution& solution) -> MineralSolutionStateFunction
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
