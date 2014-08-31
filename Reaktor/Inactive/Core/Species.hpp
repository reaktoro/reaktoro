/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Core/Types.hpp>
#include <Reaktor/Common/Units.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

/**
 * Provides a computational representation of a chemical species
 *
 * The Species class is used to represent a chemical species. It is an important
 * class in the library, since it defines fundamental attributes of a general
 * chemical species such as its elemental formula, electrical charge and molar
 * mass. In addition, it provides the functionality to calculate its standard
 * chemical potential at given temperature *T* and pressure *P* points.
 *
 * @see AqueousSpecies, GaseousSpecies, MineralSpecies, Database, ChemicalSystem
 * @ingroup Core
 */
class Species
{
public:
    /// Construct a default chemical species
    Species();

    /// Set the name of the chemical species
    auto setName(const std::string& name) -> void;

    /// Set the chemical formula of the chemical species
    auto setFormula(const std::string& formula) -> void;

    /**
     * Set the elements that constitute the chemical species
     *
     * The elements of a chemical species can be set with a sequence
     * of pairs (_element_, _number of atoms_). See below for an example
     * of setting the elements of species @f$\mathrm{CaCO_{3}(aq)}@f$:
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Species species;
     * species.setName("CaCO3(aq)");
     * species.setElements({{"Ca", 1}, {"C", 1}, {"O", 3}});
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param elements The elements and their number of atoms
     */
    auto setElements(const std::map<std::string, double>& elements) -> void;

    /**
     * Set the elements that constitute the chemical species
     *
     * The elements of a chemical species can be set with a formated string
     * @c "E1(N1)E2(N2)...Ei(Ni)", where @c Ei represents an element
     * (e.g., Na, C, H, Ca) and @c Ni the number of atoms of element @c Ei.
     * See below for an example of setting the elements of species
     * @f$\mathrm{CaCO_{3}(aq)}@f$:
     *
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Species species;
     * species.setName("CaCO3(aq)");
     * species.setElements("Ca(1)C(1)O(3)");
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * @param elements The elements and their number of atoms as a formated string
     */
    auto setElements(const std::string& elements) -> void;

    /// Set the molar mass of the chemical species
    auto setMolarMass(units::MolarMass value) -> void;

    /// Set the electrical charge of the chemical species
    auto setCharge(double value) -> void;

    /// Set the type of the chemical species
    auto setType(const std::string& type) -> void;

    /// Sets the standard chemical potential function of the chemical species
    auto setChemicalPotential(const ChemicalPotentialFn& function) -> void;

    /// Sets the standard chemical potential function of the chemical species
    auto setChemicalPotential(double constant) -> void;

    /**
     * Sets the density function of the chemical species
     * @param density The density function
     * @see DensityFn
     */
    auto setDensity(const DensityFn& density) -> void;

    /**
     * Sets the density function of the chemical species
     * @param density The constant density value
     * @see DensityFn
     */
    auto setDensity(double density) -> void;

    /// Get the name of the chemical species
    auto name() const -> const std::string&;

    /// Get the chemical formula of the species
    auto formula() const -> std::string;

    /// Get the elemental formula of the chemical species
    auto elements() const -> const std::map<std::string, double>&;

    /// Get the molar mass of the chemical species
    auto molarMass() const -> units::MolarMass;

    /// Get the electrical charge of the chemical species
    auto charge() const -> double;

    /// Get the type of the chemical species
    auto type() const -> const std::string&;

    /// Checks if the species contains an element
    auto containsElement(const std::string& element) const -> bool;

    /// Gets the number of atoms of an element in the species
    auto elementAtoms(const std::string& element) const -> double;

    /// Gets the number of atoms of the elements that compose the species
    auto elementAtoms() const -> std::vector<double>;

    /// Gets the names of the elements that compose the species
    auto elementNames() const -> std::vector<std::string>;

    /// Get the interpolation function of the standard chemical potential of the species
    auto chemicalPotential() const -> const ChemicalPotentialFn&;

    /// Calculate the standard chemical potential of the species (in units of J/mol)
    auto chemicalPotential(double T, double P) const -> double;

    /**
     * Calculates the density of the species (in units of kg/m3)
     * @param T The temperature for the calculation (in units of K)
     * @param P The pressure for the calculation (in units of bar)
     */
    auto density(double T, double P) const -> double;

    /**
     * Calculates the molar density of the species (in units of mol/m3)
     * @param T The temperature for the calculation (in units of K)
     * @param P The pressure for the calculation (in units of bar)
     */
    auto molarDensity(double T, double P) const -> double;

    /**
     * Calculates the molar volume of the species (in units of m3/mol)
     * @param T The temperature for the calculation (in units of K)
     * @param P The pressure for the calculation (in units of bar)
     */
    auto molarVolume(double T, double P) const -> double;

    /**
     * Calculates the specific volume of the species (in units of m3/kg)
     * @param T The temperature for the calculation (in units of K)
     * @param P The pressure for the calculation (in units of bar)
     */
    auto specificVolume(double T, double P) const -> double;

    /// Convert the chemical species object to a string
    operator std::string() const;

private:
    /// The name of the species
    std::string name$;

    /// The chemical formula of the species
    std::string formula$;

    /// The elemental formula of the species
    std::map<std::string, double> elements$;

    /// The molar mass of the species
    units::MolarMass molar_mass$;

    /// The electrical charge of the species
    double charge$;

    /// The type of the species (e.g., aqueous, gaseous, mineral, etc.)
    std::string type$;

    /// The standard chemical potential function of the chemical species
    ChemicalPotentialFn chemical_potential$;

    /// The density function of the chemical species
    DensityFn density$;
};

/// Convert a vector of chemical species to a vector of species names
auto names(const std::vector<Species>& species) -> std::vector<std::string>;

/// Convert a vector of chemical species to a vector of species charges
auto charges(const std::vector<Species>& species) -> Vector;

/**
 * Outputs a Species instance
 */
auto operator<<(std::ostream& out, const Species& species) -> std::ostream&;

} // namespace Reaktor
