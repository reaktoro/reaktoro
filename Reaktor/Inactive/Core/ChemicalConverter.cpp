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

#include "ChemicalConverter.hpp"


// Reaktor includes
#include <Reaktor/Core/Database.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Math/BilinearInterpolator.hpp>
#include <Reaktor/Phases/AqueousPhase.hpp>
#include <Reaktor/Phases/GaseousPhase.hpp>
#include <Reaktor/Phases/MineralPhase.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>

namespace Reaktor {
namespace internal {

auto nameMineralPhase(const MineralPhase& phase) -> std::string
{
    std::stringstream name;
    for(const MineralSpecies& iter : phase.species())
        name << iter.name() << "-";
    std::string str = name.str();
    str = str.substr(0, str.size() - 1);
    return str;
}

} /* namespace internal */

ChemicalConverter::ChemicalConverter(const Database& database)
: thermo_editor$(database)
{

}

void
ChemicalConverter::setTemperatures(const std::vector<double>& temperatures)
{
    thermo_editor$.setTemperatures(temperatures);
}

void
ChemicalConverter::setPressures(const std::vector<double>& pressures)
{
    thermo_editor$.setPressures(pressures);
}

const Database&
ChemicalConverter::database() const
{
    return thermo_editor$.database();
}

const std::vector<double>&
ChemicalConverter::temperatures() const
{
    return thermo_editor$.temperatures();
}

const std::vector<double>&
ChemicalConverter::pressures() const
{
    return thermo_editor$.pressures();
}

Species
ChemicalConverter::convertSpecies(const std::string& name) const
{
    if(database().containsAqueousSpecies(name))
        return convertSpecies(database().aqueousSpecies(name));

    if(database().containsGaseousSpecies(name))
        return convertSpecies(database().gaseousSpecies(name));

    if(database().containsMineralSpecies(name))
        return convertSpecies(database().mineralSpecies(name));

    Exception exception;
    exception.error  << "Unable to convert " << name << " to a Species instance.";
    exception.reason << "This species is not contained in the database.";
    raise(exception);

    return Species();
}

Species
ChemicalConverter::convertSpecies(const AqueousSpecies& species) const
{
    Species converted;
    converted.setName(species.name());
    converted.setCharge(species.charge());
    converted.setElements(species.elements());
    converted.setMolarMass(species.molarMass());
    converted.setChemicalPotential(thermo_editor$.thermoData(species).gibbs);
    return converted;
}

Species
ChemicalConverter::convertSpecies(const GaseousSpecies& species) const
{
    Species converted;
    converted.setName(species.name());
    converted.setCharge(0.0);
    converted.setElements(species.elements());
    converted.setMolarMass(species.molarMass());
    converted.setChemicalPotential(thermo_editor$.thermoData(species).gibbs);
    return converted;
}

Species
ChemicalConverter::convertSpecies(const MineralSpecies& species) const
{
    const double density = species.density().in(unit(kg)/unit(m3));

    Species converted;
    converted.setName(species.name());
    converted.setCharge(0.0);
    converted.setDensity(density);
    converted.setElements(species.elements());
    converted.setMolarMass(species.molarMass());
    converted.setChemicalPotential(thermo_editor$.thermoData(species).gibbs);
    return converted;
}

Phase
ChemicalConverter::convertPhase(const AqueousPhase& phase) const
{
    // Convert the aqueous species to Species instances
    std::vector<Species> species;
    for(const AqueousSpecies& iter : phase.species())
        species.push_back(convertSpecies(iter));

    // Define the concentration function of the aqueous phase
    ConcentrationFn concentration = [=](const Vector& n) -> Vector
    {
        return phase.concentrations(n);
    };

    // Define the activity function of the aqueous phase
    ActivityFn activity = [=](double T, double P, const Vector& n)
    {
        return phase.activities(T, P, n);
    };

    Phase converted;
    converted.setName("Aqueous");
    converted.setSpecies(species);
    converted.setConcentrationFn(concentration);
    converted.setActivityFn(activity);

    return converted;
}

Phase
ChemicalConverter::convertPhase(const GaseousPhase& phase) const
{
    // Convert the gaseous species to Species instances
    std::vector<Species> species;
    for(const GaseousSpecies& iter : phase.species())
        species.push_back(convertSpecies(iter));

    // Define the concentration function of the gaseous phase
    ConcentrationFn concentration = [=](const Vector& n) -> Vector
    {
        return phase.concentrations(n);
    };

    // Define the activity function of the gaseous phase
    ActivityFn activity = [=](double T, double P, const Vector& n)
    {
        return phase.activities(T, P, n);
    };

    Phase converted;
    converted.setName("Gaseous");
    converted.setSpecies(species);
    converted.setConcentrationFn(concentration);
    converted.setActivityFn(activity);

    return converted;
}

Phase
ChemicalConverter::convertPhase(const MineralPhase& phase) const
{
    // Convert the mineral species to Species instances
    std::vector<Species> species;
    for(const MineralSpecies& iter : phase.species())
        species.push_back(convertSpecies(iter));

    // Define the concentration function of the mineral phase
    ConcentrationFn concentration = [=](const Vector& n) -> Vector
    {
        return phase.concentrations(n);
    };

    // Define the activity function of the mineral phase
    ActivityFn activity = [=](double T, double P, const Vector& n)
    {
        return phase.activities(T, P, n);
    };

    Phase converted;
    converted.setName(internal::nameMineralPhase(phase));
    converted.setSpecies(species);
    converted.setConcentrationFn(concentration);
    converted.setActivityFn(activity);

    return converted;
}

} // namespace Reaktor
