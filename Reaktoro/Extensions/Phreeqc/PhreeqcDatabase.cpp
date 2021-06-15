// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "PhreeqcDatabase.hpp"

// CMakeRC includes
#include <cmrc/cmrc.hpp>

CMRC_DECLARE(ReaktoroDatabases);

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Memoization.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcThermo.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcUtils.hpp>
#include <Reaktoro/Models/ReactionThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelVantHoff.hpp>

namespace Reaktoro {
namespace detail {

/// An auxiliary type to create the Species objects from a PHREEQC database.
struct PhreeqcDatabaseHelper
{
    /// The PHREEQC instance.
    /// Note: Shared pointer ensures that the pointers in PHREEQC to
    /// species/phase objects remain valid throughout.
    SharedPtr<PHREEQC> phreeqc;

    /// The list of elements in the database
    ElementList elements;

    /// The list of species in the database
    SpeciesList species_list;

    /// The symbols of the elements already in the database (needed for fast existence check)
    Set<String> element_symbols;

    /// The names of the species already in the database (needed for fast existence check)
    Set<String> species_names;

    /// Construct a default PhreeqcDatabaseHelper object.
    PhreeqcDatabaseHelper()
    : phreeqc(new PHREEQC())
    {}

    /// Construct a PhreeqcDatabaseHelper object with given database.
    PhreeqcDatabaseHelper(String database)
    : PhreeqcDatabaseHelper()
    {
        // Load the PHREEQC database
        PhreeqcUtils::load(*phreeqc, database);

        // Create the Element objects
        for(auto i = 0; i < phreeqc->count_elements; ++i)
            addElement(phreeqc->elements[i]);

        // Create the Species objects using PhreeqcSpecies pointers (aqueous, exchange, surface species)
        for(auto i = 0; i < phreeqc->count_s; ++i)
            addSpecies(phreeqc->s[i]);

        // Create the Species objects using PhreeqcPhase pointers (gases and minerals)
        for(auto i = 0; i < phreeqc->count_phases; ++i)
            addSpecies(phreeqc->phases[i]);
    }

    /// Return true if an Element object with given symbol already exists in the database.
    auto containsElement(String symbol) -> bool
    {
        return element_symbols.find(symbol) != element_symbols.end();
    }

    /// Return true if a Species object with given name already exists in the database.
    auto containsSpecies(String name) -> bool
    {
        return species_names.find(name) != species_names.end();
    }

    /// Append an Element object with given PhreeqcElement pointer.
    auto addElement(const PhreeqcElement* e) -> void
    {
        // Skip if element with same symbol has already been appended!
        if(containsElement(PhreeqcUtils::symbol(e)))
            return;

        elements.append(Element()
            .withSymbol(PhreeqcUtils::symbol(e)) // element symbol is H, H(0), C, Fe, Fe(3), etc.
            .withName(PhreeqcUtils::name(e))     // element name is H+, H2, CO3-2, Fe+2, Fe+3, etc.
            .withMolarMass(PhreeqcUtils::molarMass(e)));

        element_symbols.insert(PhreeqcUtils::symbol(e));
    }

    /// Create a Species object with given Phreeqc species pointer.
    template<typename SpeciesType>
    auto addSpecies(const SpeciesType* s) -> void
    {
        // Skip if species with same name has already been appended!
        if(containsSpecies(PhreeqcUtils::name(s)))
            return;

        species_list.append(Species()
            .withName(PhreeqcUtils::name(s))
            .withFormula(PhreeqcUtils::formula(s))
            .withElements(createElements(s))
            .withCharge(PhreeqcUtils::charge(s))
            .withAggregateState(PhreeqcUtils::aggregateState(s))
            .withFormationReaction(createFormationReaction(s))
            .withAttachedData(s));

        species_names.insert(PhreeqcUtils::name(s));
    }

    /// Create the elements of a Species object with given Phreeqc species pointer
    template<typename SpeciesType>
    auto createElements(const SpeciesType* s) -> Pairs<Element, double>
    {
        Pairs<Element, double> pairs;
        for(auto&& [element, coeff] : PhreeqcUtils::elements(s))
        {
            const auto idx = elements.find(PhreeqcUtils::symbol(element));
            if(idx == elements.size())
                addElement(element); // create and append an Element object for this PHREEQC element first!
            pairs.emplace_back(elements[idx], coeff);
        }

        return pairs;
    }

    /// Create the formation reaction of a product species.
    template<typename SpeciesType>
    auto createFormationReaction(const SpeciesType* s) -> FormationReaction
    {
        return FormationReaction()
            .withProduct(PhreeqcUtils::name(s))
            .withReactants(createReactants(s))
            .withProductStandardVolumeModel(PhreeqcUtils::standardVolumeModel(s))
            .withReactionThermoModel(PhreeqcUtils::reactionThermoModel(s));
    }

    /// Create the reactant species and their stoichiometries in a formation reaction.
    template<typename SpeciesType>
    auto createReactants(const SpeciesType* s) -> Pairs<Species, double>
    {
        // Note: This method may introduce recursion on purpose with the
        // addSpecies method call below, when a Species object for one of its
        // reactants in its formation reaction does not yet exist.

        Pairs<Species, double> pairs;
        for(const auto& [reactant, coeff] : PhreeqcUtils::reactants(s))
        {
            const auto idx = species_list.find(PhreeqcUtils::name(reactant));
            if(idx == species_list.size())
                addSpecies(reactant); // create and append a Species object for this PHREEQC reactant species first!
            pairs.emplace_back(species_list[idx], coeff);
        }
        return pairs;
    }
};

/// Return the contents of the embedded PHREEQC database with given name (or empty)
auto getPhreeqcDatabaseContent(String name) -> String
{
    error(!oneof(name,
        "Amm.dat",
        "frezchem.dat",
        "iso.dat",
        "llnl.dat",
        "minteq.dat",
        "minteq.v4.dat",
        "phreeqc.dat",
        "pitzer.dat",
        "sit.dat",
        "wateq4f.dat"),
        "Could not load embedded PHREEQC database file with name `", name, "`. ",
        "The currently supported names are: \n"
        "    - Amm.dat       \n",
        "    - frezchem.dat  \n",
        "    - iso.dat       \n",
        "    - llnl.dat      \n",
        "    - minteq.dat    \n",
        "    - minteq.v4.dat \n",
        "    - phreeqc.dat   \n",
        "    - pitzer.dat    \n",
        "    - sit.dat       \n",
        "    - wateq4f.dat   \n",
        "");
    auto fs = cmrc::ReaktoroDatabases::get_filesystem();
    auto contents = fs.open("databases/phreeqc/" + name);
    return String(contents.begin(), contents.end());
}

/// Create the Species objects from given PHREEQC database.
auto createSpeciesWithDatabaseContentOrPath(String database)
{
    PhreeqcDatabaseHelper helper(database);
    return helper.species_list;
}

} // namespace detail

PhreeqcDatabase::PhreeqcDatabase()
: Database()
{}

PhreeqcDatabase::PhreeqcDatabase(String name)
: PhreeqcDatabase(PhreeqcDatabase::withName(name))
{}

auto PhreeqcDatabase::withName(std::string name) -> PhreeqcDatabase
{
    // PhreeqcDatabase db;
    // const auto content = detail::getPhreeqcDatabaseContent(name);
    // const auto species = detail::createSpeciesWithDatabaseContentOrPath(content);
    // db.addSpecies(species);
    // return db;
    PhreeqcDatabase db;
    const auto content = detail::getPhreeqcDatabaseContent(name);
    detail::PhreeqcDatabaseHelper helper(content);
    db.addSpecies(helper.species_list);
    db.attachData(helper);
    return db;
}

auto PhreeqcDatabase::fromFile(std::string path) -> PhreeqcDatabase
{
    // PhreeqcDatabase db;
    // const auto species = detail::createSpeciesWithDatabaseContentOrPath(path);
    // db.addSpecies(species);
    // return db;

    detail::PhreeqcDatabaseHelper helper(path);
    PhreeqcDatabase db;
    db.addSpecies(helper.species_list);
    db.attachData(helper);
    return db;
}

auto PhreeqcDatabase::load(String filename) -> PhreeqcDatabase&
{
    // Database::clear();
    // const auto species = detail::createSpeciesWithDatabaseContentOrPath(filename);
    // Database::addSpecies(species);
    // return *this;
    detail::PhreeqcDatabaseHelper helper(filename);
    Database::clear();
    Database::addSpecies(helper.species_list);
    Database::attachData(helper);
    return *this;
}

auto PhreeqcDatabase::contents(String database) -> String
{
    return detail::getPhreeqcDatabaseContent(database);
}

} // namespace Reaktoro
