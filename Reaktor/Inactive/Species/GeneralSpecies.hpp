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
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Units.hpp>

namespace Reaktor {

/**
 * Defines a base class for all species classes
 *
 * The GeneralSpecies class is used to represent the common properties of
 * chemical species. It is used as a base class for all other species classes.
 *
 * @see AqueousSpecies, GaseousSpecies, MineralSpecies, Database
 * @ingroup Species
 */
class GeneralSpecies
{
public:
    /**
     * Constructs a default species
     */
    GeneralSpecies();

    /**
     * Sets the name of the species
     */
    auto setName(const std::string& name) -> void;

    /**
     * Sets the chemical formula of the species
     */
    auto setFormula(const std::string& formula) -> void;

    /**
     * Sets the elements that constitute the species
     *
     * The elements of a species can be set with a sequence
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
     * Sets the elements that constitute the species
     *
     * The elements of a species can be set with a formated string
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

    /**
     * Sets the molar mass of the species
     *
     * Note that the molar mass of the species is set
     * by default when the elements of the species and their
     * stoichiometries are set using @ref setElements.
     *
     * @see setElements
     */
    auto setMolarMass(units::MolarMass value) -> void;

    /**
     * Sets the electrical charge of the species
     */
    auto setCharge(double value) -> void;

    /**
     * Gets the name of the species
     */
    auto name() const -> const std::string&;

    /**
     * Gets the chemical formula of the species
     */
    auto formula() const -> std::string;

    /**
     * Gets the elemental formula of the species
     */
    auto elements() const -> const std::map<std::string, double>&;

    /**
     * Gets the molar mass of the species
     */
    auto molarMass() const -> units::MolarMass;

    /**
     * Gets the electrical charge of the species
     */
    auto charge() const -> double;

    /**
     * Checks if the species contains an element
     */
    auto containsElement(const std::string& element) const -> bool;

    /**
     * Gets the number of atoms of an element in the species
     */
    auto elementAtoms(const std::string& element) const -> double;

    /**
     * Gets the number of atoms of the elements that compose the species
     */
    auto elementAtoms() const -> std::vector<double>;

    /**
     * Gets the names of the elements that compose the species
     */
    auto elementNames() const -> std::vector<std::string>;

    /**
     * Converts the species object to a string
     */
    operator std::string() const;

private:
    /// The name of the species
    std::string name$;

    /// The chemical formula of the species
    std::string formula$;

    /// The elemental formula of the species
    std::map<std::string, double> elements$;

    /// The molar mass of the species (in units of g/mol)
    units::MolarMass molar_mass$;

    /// The electrical charge of the species
    double charge$;
};

/**
 * Outputs a GeneralSpecies instance
 */
auto operator<<(std::ostream& out, const GeneralSpecies& species) -> std::ostream&;

} /* namespace Reaktor */
