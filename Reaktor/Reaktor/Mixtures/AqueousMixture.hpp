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
#include <set>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/ScalarResult.hpp>
#include <Reaktor/Common/VectorResult.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Mixtures/GeneralMixture.hpp>

namespace Reaktor {

/// Provide a computational representation of an aqueous mixture
///
/// The AqueousMixture class is defined as a collection of AqueousSpecies objects,
/// representing, therefore, a mixture of aqueous species. Its main purpose is to
/// provide the necessary operations in the calculation of activities of aqueous
/// species. It implements methods for the calculation of molar fractions, molalities,
/// stoichiometric molalities, and effective and stoichiometric ionic strengths.
/// In addition, it provides methods that retrives information about the ionic, neutral and
/// complex species.
///
/// @see AqueousSpecies
///
/// @ingroup Mixtures
class AqueousMixture : public GeneralMixture<AqueousSpecies>
{
public:
    /// Construct a default AqueousMixture instance
    AqueousMixture();

    /// Construct an AqueousMixture instance with given species
    /// @param species The species that compose the aqueous mixture
    AqueousMixture(const std::vector<AqueousSpecies>& species);

    /// Destroy the AqueousMixture instance
    virtual ~AqueousMixture();

    /// Get the index of the water species H<sub>2</sub>O(l)
    auto indexWater() const -> Index;

    /// Get the indices of the charged aqueous species in the aqueous mixture
    auto indicesChargedSpecies() const -> const Indices&;

    /// Get the indices of the neutral aqueous species in the aqueous mixture
    auto indicesNeutralSpecies() const -> const Indices&;

    /// Get the electrical charges of the aqueous species
	auto charges() const -> Vector;

	/// Get the electrical charges of the charged species only
	auto chargesChargedSpecies() const -> Vector;

    /// Get the dissociation matrix of the neutral species into charged species
    ///
    /// This matrix defines the stoichiometric relationship between the neutral
    /// and charged species produced from their dissociation. For example, the
    /// stoichiometry of the j-th charged species in the dissociation of the i-th
    /// neutral species is given by the (i, j)-th entry in the dissociation matrix.
    auto dissociationMatrix() const -> const Matrix&;

private:
    /// The index of the water species
    Index index_water;

    /// The indices of the neutral aqueous species
    Indices indices_neutral_species;

    /// The indices of the charged aqueous species
    Indices indices_charged_species;

    /// The electrical charges of the aqueous species
    Vector z;

    /// The electrical charges of the ionic components
    Vector z_charged;

    /// The matrix that represents the dissociation of the neutral species into charged species
    Matrix dissociation_matrix;
};

/// Get the number of charged aqueous species in an aqueous mixture
/// @param mixture The aqueous mixture
auto numChargedSpecies(const AqueousMixture& mixture) -> unsigned;

/// Get the number of ions in an aqueous mixture
/// @param mixture The aqueous mixture
auto numNeutralSpecies(const AqueousMixture& mixture) -> unsigned;

/// Get the local index of a given charged species among the charged species in an aqueous mixture
/// @param mixture The aqueous mixture
/// @param name The name of the charged species
/// @return The local index of the charged species if found, or the number of charged species otherwise
auto localIndexChargedSpecies(const AqueousMixture& mixture, const std::string& name) -> Index;

/// Get the local index of a given neutral species among the neutral species in an aqueous mixture
/// @param mixture The aqueous mixture
/// @param name The name of the neutral species
/// @return The local index of the neutral species if found, or the number of neutral species otherwise
auto localIndexNeutralSpecies(const AqueousMixture& mixture, const std::string& name) -> Index;

/// Get the local index of a given cation among the cations in an aqueous mixture
/// @param mixture The aqueous mixture
/// @param name The name of the cation
/// @return The local index of the cation if found, or the number of cations otherwise
auto localIndexCation(const AqueousMixture& mixture, const std::string& name) -> Index;

/// Get the local index of a given anion among the anions in an aqueous mixture
/// @param mixture The aqueous mixture
/// @param name The name of the anion
/// @return The local index of the anion if found, or the number of anions otherwise
auto localIndexAnion(const AqueousMixture& mixture, const std::string& name) -> Index;

/// Get the indices of the cations in an aqueous mixture
/// @param mixture The aqueous mixture
auto indicesCations(const AqueousMixture& mixture) -> Indices;

/// Get the indices of the anions in an aqueous mixture
/// @param mixture The aqueous mixture
auto indicesAnions(const AqueousMixture& mixture) -> Indices;

/// Get the names of the charged species in an aqueous mixture
/// @param mixture The aqueous mixture
auto namesChargedSpecies(const AqueousMixture& mixture) -> std::vector<std::string>;

/// Get the names of the neutral species in an aqueous mixture
/// @param mixture The aqueous mixture
auto namesNeutralSpecies(const AqueousMixture& mixture) -> std::vector<std::string>;

/// Get the names of the cations in an aqueous mixture
/// @param mixture The aqueous mixture
auto namesCations(const AqueousMixture& mixture) -> std::vector<std::string>;

/// Get the names of the cations in an aqueous mixture
/// @param mixture The aqueous mixture
auto namesAnions(const AqueousMixture& mixture) -> std::vector<std::string>;

/// Calculate the molalities of the aqueous species and its molar derivatives
/// @param mixture The aqueous mixture
/// @param n The molar abundance of species (in units of mol)
/// @return The molalities and their molar derivatives
auto molalities(const AqueousMixture& mixture, const Vector& n) -> VectorResult;

/// Calculate the stoichiometric molalities of the ions and its molar derivatives
/// @param mixture The aqueous mixture
/// @param m The molalities of the aqueous species and their molar derivatives
/// @return The stoichiometric molalities and their molar derivatives
auto molalitiesStoichiometric(const AqueousMixture& mixture, const VectorResult& m) -> VectorResult;

/// Calculate the effective ionic strength of the aqueous mixture and its molar derivatives
/// @param mixture The aqueous mixture
/// @param m The molalities of the aqueous species and their molar derivatives
/// @return The effective ionic strength of the aqueous mixture and its molar derivatives
auto ionicStrength(const AqueousMixture& mixture, const VectorResult& m) -> ScalarResult;

/// Calculate the stoichiometric ionic strength of the aqueous mixture and its molar derivatives
/// @param mixture The aqueous mixture
/// @param ms The stoichiometric molalities of the ions and their molar derivatives
/// @return The stoichiometric ionic strength of the aqueous mixture and its molar derivatives
auto ionicStrengthStoichiometric(const AqueousMixture& mixture, const VectorResult& ms) -> ScalarResult;

} // namespace Reaktor
