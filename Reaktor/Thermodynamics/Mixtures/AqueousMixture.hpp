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
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktor/Thermodynamics/Mixtures/GeneralMixture.hpp>

namespace Reaktor {

/**
 * Provides a computational representation of an aqueous mixture
 *
 * The AqueousMixture class is defined as a collection of AqueousSpecies objects,
 * representing, therefore, a mixture of aqueous species. Its main purpose is to
 * provide the necessary operations in the calculation of activities of aqueous
 * species. It implements methods for the calculation of molar fractions, molalities,
 * stoichiometric molalities, and effective and stoichiometric ionic strengths.
 * In addition, it provides methods that retrives information about the ionic, neutral and
 * complex species.
 *
 * @see AqueousSpecies
 *
 * @ingroup Mixtures
 */
class AqueousMixture : public GeneralMixture<AqueousSpecies>
{
public:
    /**
     * Constructs a default AqueousMixture instance
     */
    AqueousMixture();

    /**
     * Constructs an AqueousMixture instance with given species
     * @param species The species that compose the aqueous mixture
     */
    AqueousMixture(const std::vector<AqueousSpecies>& species);

    /**
     * Destroys the instance
     */
    virtual ~AqueousMixture();

    /**
     * Gets the number of charged aqueous species in the aqueous mixture
     */
    auto numChargedSpecies() const -> unsigned;

    /**
     * Gets the number of ions in the aqueous mixture
     */
    auto numIons() const -> unsigned;

    /**
     * Gets the number of aqueous complexes in the aqueous mixture
     */
    auto numComplexes() const -> unsigned;

    /**
     * Gets the electrical charges of the aqueous species
     */
    auto charges() const -> Vector;

    /**
     * Gets the indices of the neutral aqueous species in the aqueous mixture
     */
    auto idxNeutralSpecies() const -> const Indices&;

    /**
     * Gets the indices of the charged aqueous species in the aqueous mixture
     *
     * The charged aqueous species are defined as aqueous species that possess
     * non-zero electrical charges.
     */
    auto idxChargedSpecies() const -> const Indices&;

    /**
     * Gets the indices of the ions in the aqueous mixture
     *
     * The ions are defined as charged aqueous species that are produced from
     * the dissociation of aqueous complexes. For example, if the dissociation
     * of the species @f$\ce{AgCl2-}@f$ into @f$\ce{Ag+}@f$ and @f$\ce{Cl-}@f$
     * is provided, then @f$\ce{AgCl2-}@f$ will be considered a charged aqueous
     * species, and @f$\ce{Ag+}@f$ and @f$\ce{Cl-}@f$ ionic species.
     */
    auto idxIons() const -> const Indices&;

    /**
     * Gets the indices of the cations in the aqueous mixture
     */
    auto idxCations() const -> Indices;

    /**
     * Gets the indices of the anions in the aqueous mixture
     */
    auto idxAnions() const -> Indices;

    /**
     * Gets the indices of the aqueous complexes in the aqueous mixture
     *
     * The aqueous complexes are definition as aqueous species that can dissociate into elementary
     * ions. These can be charged or neutral, such as @f$\ce{NaCl(aq)}@f$, @f$\ce{HCl(aq)}@f$ and
     * @f$\ce{CaOH+}@f$.
     */
    auto idxComplexes() const -> const Indices&;

    /**
     * Gets the index of the water species @f$\ce{H2O(l)}@f$.
     */
    auto idxWater() const -> const Index&;

    /**
     * Gets the local index of an ion among the ions in the aqueous mixture
     * @param ion The name of the ion to be found among the ions in the mixture
     * @return The local index of the ion if found. The number of ions otherwise.
     */
    auto idxIon(const std::string& ion) const -> Index;

    /**
     * Gets the names of the neutral species in the aqueous mixture
     */
    auto neutralSpecies() const -> std::vector<std::string>;

    /**
     * Gets the names of the charged species in the aqueous mixture
     */
    auto chargedSpecies() const -> std::vector<std::string>;

    /**
     * Gets the names of the cations in the aqueous mixture
     */
    auto cations() const -> std::vector<std::string>;

    /**
     * Gets the names of the cations in the aqueous mixture
     */
    auto anions() const -> std::vector<std::string>;

    /**
     * Gets the names of the complex species in the aqueous mixture
     */
    auto complexes() const -> std::vector<std::string>;

    /**
     * Gets the dissociation matrix of the aqueous complexes into ions
     *
     * This the matrix defines the stoichiometric relationship between the aqueous complexes and the
     * ions produced from their dissociation. For example, the stoichiometry of the *j*-th ion in
     * the dissociation reaction of the *i*-th aqueous complex is given by the (*i*, *j*)-th entry in
     * the matrix.
     */
    auto dissociationMatrix() const -> const Matrix&;

    /**
     * Calculates the molalities of the aqueous species and its molar derivatives
     * @param n The molar abundance of species (in units of mol)
     * @return The molalities and their molar derivatives
     */
    auto molalities(const Vector& n) const -> ChemicalVector;

    /**
     * Calculates the stoichiometric molalities of the ions and its molar derivatives
     * @param m The molalities of the aqueous species and their molar derivatives
     * @return The stoichiometric molalities and their molar derivatives
     */
    auto stoichiometricMolalities(const ChemicalVector& m) const -> ChemicalVector;

    /**
     * Calculates the effective ionic strength of the aqueous mixture and its molar derivatives
     * @param m The molalities of the aqueous species and their molar derivatives
     * @return The effective ionic strength of the aqueous mixture and its molar derivatives
     */
    auto effectiveIonicStrength(const ChemicalVector& m) const -> ChemicalScalar;

    /**
     * Calculates the stoichiometric ionic strength of the aqueous mixture and its molar derivatives
     * @param ms The stoichiometric molalities of the ions and their molar derivatives
     * @return The stoichiometric ionic strength of the aqueous mixture and its molar derivatives
     */
    auto stoichiometricIonicStrength(const ChemicalVector& ms) const -> ChemicalScalar;

private:
    /// The index of the water species
    Index idx_water;

    /// The indices of the neutral aqueous species
    Indices idx_neutral_species;

    /// The indices of the charged aqueous species
    Indices idx_charged;

    /// The indices of the ionic components in the aqueous mixture.
    Indices idx_ions;

    /// The indices of the aqueous complexes in the aqueous mixture
    Indices idx_complexes;

    /// The electrical charges of the aqueous species
    Vector z;

    /// The electrical charges of the ionic components
    Vector zi;

    /// The matrix that represents the dissociation of the aqueous complexes into ions
    Matrix nu;

    /// The set of aqueous species that are regarded as ions in the aqueous mixture
    std::set<std::string> ions;
};

} // namespace Reaktor
