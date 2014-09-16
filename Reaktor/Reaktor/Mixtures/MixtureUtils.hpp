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
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>

namespace Reaktor {

typedef std::vector<AqueousSpecies> AqueousMixture;
typedef std::vector<GaseousSpecies> GaseousMixture;
typedef std::vector<MineralSpecies> MineralMixture;

/// Get the number of species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
template<typename Mixture>
auto numSpecies(const Mixture& mixture) -> unsigned;

/// Get the index of a species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
/// @param name The name of the species in the mixture
/// @return The index of the species if found, or the number of species otherwise
template<class Mixture>
auto indexSpecies(const Mixture& mixture, const std::string& name) -> Index;

/// Get the index of the water species H<sub>2</sub>O(l) in an aqueous mixture
/// @param mixture The aqueous mixture
auto indexWater(const AqueousMixture& mixture) -> Index;

/// Get the indices of the charged aqueous species in the aqueous mixture
/// @param mixture The aqueous mixture
auto indicesChargedSpecies(const AqueousMixture& mixture) -> Indices;

/// Get the indices of the neutral aqueous species in the aqueous mixture
/// @param mixture The aqueous mixture
auto indicesNeutralSpecies(const AqueousMixture& mixture) -> Indices;

/// Get the indices of the cations in an aqueous mixture
/// @param mixture The aqueous mixture
auto indicesCations(const AqueousMixture& mixture) -> Indices;

/// Get the indices of the anions in an aqueous mixture
/// @param mixture The aqueous mixture
auto indicesAnions(const AqueousMixture& mixture) -> Indices;

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

/// Get the names of the species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
template<typename Mixture>
auto namesSpecies(const Mixture& mixture) -> std::vector<std::string>;

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

/// Get the charges of the species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
template<typename Mixture>
auto chargesSpecies(const Mixture& mixture) -> Vector;

/// Get the electrical charges of the charged species
auto chargesChargedSpecies(const AqueousMixture& mixture) -> Vector;

/// Get the electrical charges of the cations
auto chargesCations(const AqueousMixture& mixture) -> Vector;

/// Get the electrical charges of the anions
auto chargesAnions(const AqueousMixture& mixture) -> Vector;

/// Get the dissociation matrix of the neutral species into charged species.
/// This matrix defines the stoichiometric relationship between the neutral
/// and charged species produced from their dissociation. For example, the
/// stoichiometry of the j-th charged species in the dissociation of the i-th
/// neutral species is given by the (i, j)-th entry in the dissociation matrix.
auto dissociationMatrix(const AqueousMixture& mixture) -> Matrix;

/// Calculate the molar fractions of the species and its molar derivatives
/// @param n The molar abundance of the species (in units of mol)
/// @return The molar fractions and its molar derivatives
auto molarFractions(const Vector& n) -> ThermoVector;

/// Calculate the molalities of the aqueous species and its molar derivatives
/// @param mixture The aqueous mixture
/// @param n The molar abundance of species (in units of mol)
/// @return The molalities and their molar derivatives
auto molalities(const AqueousMixture& mixture, const Vector& n) -> ThermoVector;

/// Calculate the stoichiometric molalities of the ions and its molar derivatives
/// @param mixture The aqueous mixture
/// @param m The molalities of the aqueous species and their molar derivatives
/// @return The stoichiometric molalities and their molar derivatives
auto molalitiesStoichiometric(const AqueousMixture& mixture, const ThermoVector& m) -> ThermoVector;

/// Calculate the effective ionic strength of the aqueous mixture and its molar derivatives
/// @param mixture The aqueous mixture
/// @param m The molalities of the aqueous species and their molar derivatives
/// @return The effective ionic strength of the aqueous mixture and its molar derivatives
auto ionicStrength(const AqueousMixture& mixture, const ThermoVector& m) -> ThermoScalar;

/// Calculate the stoichiometric ionic strength of the aqueous mixture and its molar derivatives
/// @param mixture The aqueous mixture
/// @param ms The stoichiometric molalities of the ions and their molar derivatives
/// @return The stoichiometric ionic strength of the aqueous mixture and its molar derivatives
auto ionicStrengthStoichiometric(const AqueousMixture& mixture, const ThermoVector& ms) -> ThermoScalar;

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
    ThermoVector x;
};

/// Compare two MixtureState instances for equality
auto operator==(const MixtureState& l, const MixtureState& r) -> bool;

/// A type used to describe the state of an aqueous mixture
struct AqueousMixtureState : public MixtureState
{
    /// The effective ionic strength of the aqueous mixture and its molar derivatives (in units of mol/kg)
    ThermoScalar Ie;

    /// The stoichiometric ionic strength of the aqueous mixture and its molar derivatives (in units of mol/kg)
    ThermoScalar Is;

    /// The molalities of the aqueous species and its molar derivatives (in units of mol/kg)
    ThermoVector m;

    /// The stoichiometric molalities of the ionic species and its molar derivatives (in units of mol/kg)
    ThermoVector ms;
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

auto aqueousMixtureStateFunction(const AqueousMixture& mixture) -> AqueousMixtureStateFunction;

auto gaseousMixtureStateFunction(const GaseousMixture& mixture) -> GaseousMixtureStateFunction;

auto mineralMixtureStateFunction(const MineralMixture& mixture) -> MineralMixtureStateFunction;

} // namespace Reaktor

#include "MixtureUtils.hxx"
