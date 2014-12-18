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

#pragma once

// C++ includes
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>

namespace Reaktor {

/// A type used to describe a collection of AqueousSpecies instances
typedef std::vector<AqueousSpecies> AqueousSolution;

/// A type used to describe a collection of GaseousSpecies instances
typedef std::vector<GaseousSpecies> GaseousSolution;

/// A type used to describe a collection of MineralSpecies instances
typedef std::vector<MineralSpecies> MineralSolution;

/// A type used to describe the state of an aqueous solution
struct SolutionState
{
    /// The temperature of the aqueous solution (in units of K)
    double T;

    /// The pressure of the aqueous solution (in units of Pa)
    double P;

    /// The amounts of the species in the aqueous solution (in units of mol)
    Vector n;

    /// The molar fractions of the aqueous species and its molar derivatives
    ChemicalVector x;
};

/// Compare two SolutionState instances for equality
auto operator==(const SolutionState& l, const SolutionState& r) -> bool;

/// A type used to describe the state of an aqueous solution
struct AqueousSolutionState : public SolutionState
{
    /// The effective ionic strength of the aqueous solution and its molar derivatives (in units of mol/kg)
    ChemicalScalar Ie;

    /// The stoichiometric ionic strength of the aqueous solution and its molar derivatives (in units of mol/kg)
    ChemicalScalar Is;

    /// The molalities of the aqueous species and its molar derivatives (in units of mol/kg)
    ChemicalVector m;

    /// The stoichiometric molalities of the ionic species and its molar derivatives (in units of mol/kg)
    ChemicalVector ms;
};

/// A type used to describe the state of a gaseous solution
struct GaseousSolutionState : public SolutionState
{};

/// A type used to describe the state of a mineral solution
struct MineralSolutionState : public SolutionState
{};

/// A type that describes a function that computes the state of an aqueous solution
typedef std::function<
    AqueousSolutionState(double T, double P, const Vector& n)>
        AqueousSolutionStateFunction;

/// A type that describes a function that computes the state of a gaseous solution
typedef std::function<
    GaseousSolutionState(double T, double P, const Vector& n)>
        GaseousSolutionStateFunction;

/// A type that describes a function that computes the state of a mineral solution
typedef std::function<
    MineralSolutionState(double T, double P, const Vector& n)>
        MineralSolutionStateFunction;

/// Get the number of species in a solution
/// @param solution The solution (e.g., aqueous, gaseous, mineral, etc.)
template<typename Solution>
auto numSpecies(const Solution& solution) -> unsigned;

/// Get the index of a species in a solution
/// @param solution The solution (e.g., aqueous, gaseous, mineral, etc.)
/// @param name The name of the species in the solution
/// @return The index of the species if found, or the number of species otherwise
template<class Solution>
auto speciesIndex(const Solution& solution, const std::string& name) -> Index;

/// Get the names of the species in a solution
/// @param solution The solution (e.g., aqueous, gaseous, mineral, etc.)
template<typename Solution>
auto speciesNames(const Solution& solution) -> std::vector<std::string>;

/// Get the charges of the species in a solution
/// @param solution The solution (e.g., aqueous, gaseous, mineral, etc.)
template<typename Solution>
auto speciesCharges(const Solution& solution) -> Vector;

/// Get the indices of the charged aqueous species in the aqueous solution
/// @param solution The aqueous solution
auto chargedSpeciesIndices(const AqueousSolution& solution) -> Indices;

/// Get the local index of a given charged species among the charged species in an aqueous solution
/// @param solution The aqueous solution
/// @param name The name of the charged species
/// @return The local index of the charged species if found, or the number of charged species otherwise
auto chargedSpeciesLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index;

/// Get the names of the charged species in an aqueous solution
/// @param solution The aqueous solution
auto chargedSpeciesNames(const AqueousSolution& solution) -> std::vector<std::string>;

/// Get the electrical charges of the charged species
/// @param solution The aqueous solution
auto chargedSpeciesCharges(const AqueousSolution& solution) -> Vector;

/// Get the indices of the neutral aqueous species in the aqueous solution
/// @param solution The aqueous solution
auto neutralSpeciesIndices(const AqueousSolution& solution) -> Indices;

/// Get the local index of a given neutral species among the neutral species in an aqueous solution
/// @param solution The aqueous solution
/// @param name The name of the neutral species
/// @return The local index of the neutral species if found, or the number of neutral species otherwise
auto neutralSpeciesLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index;

/// Get the names of the neutral species in an aqueous solution
/// @param solution The aqueous solution
auto neutralSpeciesNames(const AqueousSolution& solution) -> std::vector<std::string>;

/// Get the indices of the cations in an aqueous solution
/// @param solution The aqueous solution
auto cationIndices(const AqueousSolution& solution) -> Indices;

/// Get the local index of a given cation among the cations in an aqueous solution
/// @param solution The aqueous solution
/// @param name The name of the cation
/// @return The local index of the cation if found, or the number of cations otherwise
auto cationLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index;

/// Get the names of the cations in an aqueous solution
/// @param solution The aqueous solution
auto cationNames(const AqueousSolution& solution) -> std::vector<std::string>;

/// Get the electrical charges of the cations
/// @param solution The aqueous solution
auto cationCharges(const AqueousSolution& solution) -> Vector;

/// Get the indices of the anions in an aqueous solution
/// @param solution The aqueous solution
auto anionIndices(const AqueousSolution& solution) -> Indices;

/// Get the local index of a given anion among the anions in an aqueous solution
/// @param solution The aqueous solution
/// @param name The name of the anion
/// @return The local index of the anion if found, or the number of anions otherwise
auto anionLocalIndex(const AqueousSolution& solution, const std::string& name) -> Index;

/// Get the names of the cations in an aqueous solution
/// @param solution The aqueous solution
auto anionNames(const AqueousSolution& solution) -> std::vector<std::string>;

/// Get the electrical charges of the anions
/// @param solution The aqueous solution
auto anionCharges(const AqueousSolution& solution) -> Vector;

/// Get the index of the water species H<sub>2</sub>O(l) in an aqueous solution
/// @param solution The aqueous solution
auto waterIndex(const AqueousSolution& solution) -> Index;

/// Get the dissociation matrix of the neutral species into charged species.
/// This matrix defines the stoichiometric relationship between the neutral
/// and charged species produced from their dissociation. For example, the
/// stoichiometry of the j-th charged species in the dissociation of the i-th
/// neutral species is given by the (i, j)-th entry in the dissociation matrix.
/// @param solution The aqueous solution
auto dissociationMatrix(const AqueousSolution& solution) -> Matrix;

/// Calculate the molar fractions from the molar amounts of a set of species
/// @param n The molar amounts of the species
auto molarFractions(const Vector& n) -> ChemicalVector;

/// Create a function for the state calculation of an aqueous solution
/// @param solution The aqueous solution
auto aqueousSolutionStateFunction(const AqueousSolution& solution) -> AqueousSolutionStateFunction;

/// Create a function for the state calculation of a gaseous solution
/// @param solution The gaseous solution
auto gaseousSolutionStateFunction(const GaseousSolution& solution) -> GaseousSolutionStateFunction;

/// Create a function for the state calculation of a mineral solution
/// @param solution The mineral solution
auto mineralSolutionStateFunction(const MineralSolution& solution) -> MineralSolutionStateFunction;

} // namespace Reaktor

#include "ActivityUtils.hxx"
