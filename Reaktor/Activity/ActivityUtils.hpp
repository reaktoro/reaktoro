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
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>

namespace Reaktor {

/// A type used to describe a collection of AqueousSpecies instances
typedef std::vector<AqueousSpecies> AqueousMixture;

/// A type used to describe a collection of GaseousSpecies instances
typedef std::vector<GaseousSpecies> GaseousMixture;

/// A type used to describe a collection of MineralSpecies instances
typedef std::vector<MineralSpecies> MineralMixture;

/// A type used to describe the state of an aqueous mixture
struct MixtureState
{
    /// The temperature of the aqueous mixture (in units of K)
    double T;

    /// The pressure of the aqueous mixture (in units of Pa)
    double P;

    /// The amounts of the species in the aqueous mixture (in units of mol)
    Vector n;

    /// The molar fractions of the aqueous species and its molar derivatives
    ChemicalVector x;
};

/// Compare two MixtureState instances for equality
auto operator==(const MixtureState& l, const MixtureState& r) -> bool;

/// A type used to describe the state of an aqueous mixture
struct AqueousMixtureState : public MixtureState
{
    /// The effective ionic strength of the aqueous mixture and its molar derivatives (in units of mol/kg)
    ChemicalScalar Ie;

    /// The stoichiometric ionic strength of the aqueous mixture and its molar derivatives (in units of mol/kg)
    ChemicalScalar Is;

    /// The molalities of the aqueous species and its molar derivatives (in units of mol/kg)
    ChemicalVector m;

    /// The stoichiometric molalities of the ionic species and its molar derivatives (in units of mol/kg)
    ChemicalVector ms;
};

/// A type used to describe the state of a gaseous mixture
struct GaseousMixtureState : public MixtureState
{};

/// A type used to describe the state of a mineral mixture
struct MineralMixtureState : public MixtureState
{};

/// A type that describes a function that computes the state of an aqueous mixture
typedef std::function<
    AqueousMixtureState(double T, double P, const Vector& n)>
        AqueousMixtureStateFunction;

/// A type that describes a function that computes the state of a gaseous mixture
typedef std::function<
    GaseousMixtureState(double T, double P, const Vector& n)>
        GaseousMixtureStateFunction;

/// A type that describes a function that computes the state of a mineral mixture
typedef std::function<
    MineralMixtureState(double T, double P, const Vector& n)>
        MineralMixtureStateFunction;

/// Get the number of species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
template<typename Mixture>
auto numSpecies(const Mixture& mixture) -> unsigned;

/// Get the index of a species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
/// @param name The name of the species in the mixture
/// @return The index of the species if found, or the number of species otherwise
template<class Mixture>
auto speciesIndex(const Mixture& mixture, const std::string& name) -> Index;

/// Get the names of the species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
template<typename Mixture>
auto speciesNames(const Mixture& mixture) -> std::vector<std::string>;

/// Get the charges of the species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
template<typename Mixture>
auto speciesCharges(const Mixture& mixture) -> Vector;

/// Get the indices of the charged aqueous species in the aqueous mixture
/// @param mixture The aqueous mixture
auto chargedSpeciesIndices(const AqueousMixture& mixture) -> Indices;

/// Get the local index of a given charged species among the charged species in an aqueous mixture
/// @param mixture The aqueous mixture
/// @param name The name of the charged species
/// @return The local index of the charged species if found, or the number of charged species otherwise
auto chargedSpeciesLocalIndex(const AqueousMixture& mixture, const std::string& name) -> Index;

/// Get the names of the charged species in an aqueous mixture
/// @param mixture The aqueous mixture
auto chargedSpeciesNames(const AqueousMixture& mixture) -> std::vector<std::string>;

/// Get the electrical charges of the charged species
/// @param mixture The aqueous mixture
auto chargedSpeciesCharges(const AqueousMixture& mixture) -> Vector;

/// Get the indices of the neutral aqueous species in the aqueous mixture
/// @param mixture The aqueous mixture
auto neutralSpeciesIndices(const AqueousMixture& mixture) -> Indices;

/// Get the local index of a given neutral species among the neutral species in an aqueous mixture
/// @param mixture The aqueous mixture
/// @param name The name of the neutral species
/// @return The local index of the neutral species if found, or the number of neutral species otherwise
auto neutralSpeciesLocalIndex(const AqueousMixture& mixture, const std::string& name) -> Index;

/// Get the names of the neutral species in an aqueous mixture
/// @param mixture The aqueous mixture
auto neutralSpeciesNames(const AqueousMixture& mixture) -> std::vector<std::string>;

/// Get the indices of the cations in an aqueous mixture
/// @param mixture The aqueous mixture
auto cationIndices(const AqueousMixture& mixture) -> Indices;

/// Get the local index of a given cation among the cations in an aqueous mixture
/// @param mixture The aqueous mixture
/// @param name The name of the cation
/// @return The local index of the cation if found, or the number of cations otherwise
auto cationLocalIndex(const AqueousMixture& mixture, const std::string& name) -> Index;

/// Get the names of the cations in an aqueous mixture
/// @param mixture The aqueous mixture
auto cationNames(const AqueousMixture& mixture) -> std::vector<std::string>;

/// Get the electrical charges of the cations
/// @param mixture The aqueous mixture
auto cationCharges(const AqueousMixture& mixture) -> Vector;

/// Get the indices of the anions in an aqueous mixture
/// @param mixture The aqueous mixture
auto anionIndices(const AqueousMixture& mixture) -> Indices;

/// Get the local index of a given anion among the anions in an aqueous mixture
/// @param mixture The aqueous mixture
/// @param name The name of the anion
/// @return The local index of the anion if found, or the number of anions otherwise
auto anionLocalIndex(const AqueousMixture& mixture, const std::string& name) -> Index;

/// Get the names of the cations in an aqueous mixture
/// @param mixture The aqueous mixture
auto anionNames(const AqueousMixture& mixture) -> std::vector<std::string>;

/// Get the electrical charges of the anions
/// @param mixture The aqueous mixture
auto anionCharges(const AqueousMixture& mixture) -> Vector;

/// Get the index of the water species H<sub>2</sub>O(l) in an aqueous mixture
/// @param mixture The aqueous mixture
auto waterIndex(const AqueousMixture& mixture) -> Index;

/// Get the dissociation matrix of the neutral species into charged species.
/// This matrix defines the stoichiometric relationship between the neutral
/// and charged species produced from their dissociation. For example, the
/// stoichiometry of the j-th charged species in the dissociation of the i-th
/// neutral species is given by the (i, j)-th entry in the dissociation matrix.
/// @param mixture The aqueous mixture
auto dissociationMatrix(const AqueousMixture& mixture) -> Matrix;

/// Create a function for the state calculation of an aqueous mixture
/// @param mixture The aqueous mixture
auto aqueousMixtureStateFunction(const AqueousMixture& mixture) -> AqueousMixtureStateFunction;

/// Create a function for the state calculation of a gaseous mixture
/// @param mixture The gaseous mixture
auto gaseousMixtureStateFunction(const GaseousMixture& mixture) -> GaseousMixtureStateFunction;

/// Create a function for the state calculation of a mineral mixture
/// @param mixture The mineral mixture
auto mineralMixtureStateFunction(const MineralMixture& mixture) -> MineralMixtureStateFunction;

} // namespace Reaktor

#include "ActivityUtils.hxx"
