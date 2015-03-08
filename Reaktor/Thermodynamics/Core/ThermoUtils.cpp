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

#include "ThermoUtils.hpp"

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Common/Optional.hpp>
#include <Reaktor/Thermodynamics/Core/Database.hpp>
#include <Reaktor/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesThermoState.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesThermoStateHKF.hpp>

namespace Reaktor {
namespace {

auto errorNonExistentSpecies(const std::string& name) -> void
{
    Exception exception;
    exception.error << "Cannot get an instance of the species named " << name << " in the database.";
    exception.reason << "There is no such species in the database.";
    RaiseError(exception);
}

auto getSpeciesThermoProperties(std::string species, const Database& database) -> Optional<SpeciesThermoProperties>
{
    if(database.containsAqueousSpecies(species))
        return database.aqueousSpecies(species).thermo.properties;
    if(database.containsGaseousSpecies(species))
        return database.gaseousSpecies(species).thermo.properties;
    if(database.containsMineralSpecies(species))
        return database.mineralSpecies(species).thermo.properties;
    errorNonExistentSpecies(species);
    return {};
}

auto getReactionThermoProperties(std::string species, const Database& database) -> Optional<ReactionThermoProperties>
{
    if(database.containsAqueousSpecies(species))
        return database.aqueousSpecies(species).thermo.reaction;
    if(database.containsGaseousSpecies(species))
        return database.gaseousSpecies(species).thermo.reaction;
    if(database.containsMineralSpecies(species))
        return database.mineralSpecies(species).thermo.reaction;
    errorNonExistentSpecies(species);
    return {};
}

auto hasThermoParamsHKF(std::string species, const Database& database) -> bool
{
    if(database.containsAqueousSpecies(species))
        return not database.aqueousSpecies(species).thermo.hkf.empty();
    if(database.containsGaseousSpecies(species))
        return not database.gaseousSpecies(species).thermo.hkf.empty();
    if(database.containsMineralSpecies(species))
        return not database.mineralSpecies(species).thermo.hkf.empty();
    errorNonExistentSpecies(species);
    return {};
}

auto speciesThermoStateHKF(double T, double P, std::string species, const Database& database) -> SpeciesThermoState
{
    if(database.containsAqueousSpecies(species))
        return speciesThermoStateHKF(T, P, database.aqueousSpecies(species));
    if(database.containsGaseousSpecies(species))
        return speciesThermoStateHKF(T, P, database.gaseousSpecies(species));
    if(database.containsMineralSpecies(species))
        return speciesThermoStateHKF(T, P, database.mineralSpecies(species));
    errorNonExistentSpecies(species);
    return {};
}

template<typename PropertyFunction, typename EvalFunction>
auto standardPropertyFromReaction(double T, double P, std::string species, const Database& database,
    const ReactionThermoProperties& reaction, PropertyFunction property, EvalFunction eval) -> ThermoScalar
{
    double sum = 0.0;
    for(auto pair : reaction.equation)
    {
        const auto reactant = pair.first;
        const auto stoichiometry = pair.second;
        if(reactant != species)
            sum -= stoichiometry * property(T, P, reactant, database).val();
    }
    const double stoichiometry = reaction.equation.at(species);
    sum += eval();
    sum /= stoichiometry;
    return {sum, 0.0, 0.0};
}

auto standardGibbsEnergyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.gibbs_energy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, standardGibbsEnergy, eval);
}

auto standardHelmholtzEnergyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.helmholtz_energy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, standardHelmholtzEnergy, eval);
}

auto standardInternalEnergyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.internal_energy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, standardInternalEnergy, eval);
}

auto standardEnthalpyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.enthalpy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, standardEnthalpy, eval);
}

auto standardEntropyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.entropy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, standardEntropy, eval);
}

auto standardVolumeFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.volume(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, standardVolume, eval);
}

auto standardHeatCapacityFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.heat_capacity(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, standardHeatCapacity, eval);
}

} // namespace

auto standardGibbsEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    const auto species_thermo_properties = getSpeciesThermoProperties(species, database);
    if(not species_thermo_properties.empty())
        if(not species_thermo_properties().gibbs_energy.empty())
            return ThermoScalar(species_thermo_properties().gibbs_energy(T, P), 0.0, 0.0);

    const auto reaction_thermo_properties = getReactionThermoProperties(species, database);
    if(not reaction_thermo_properties.empty())
        if(not reaction_thermo_properties().gibbs_energy.empty())
            return standardGibbsEnergyFromReaction(T, P, species, database, reaction_thermo_properties());

    if(hasThermoParamsHKF(species, database))
        return speciesThermoStateHKF(T, P, species, database).gibbs_energy;

    Exception exception;
    exception.error << "Cannot calculate the standard Gibbs energy of species " << species << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

auto standardHelmholtzEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    const auto species_thermo_properties = getSpeciesThermoProperties(species, database);
    if(not species_thermo_properties.empty())
        if(not species_thermo_properties().helmholtz_energy.empty())
            return ThermoScalar(species_thermo_properties().helmholtz_energy(T, P), 0.0, 0.0);

    const auto reaction_thermo_properties = getReactionThermoProperties(species, database);
    if(not reaction_thermo_properties.empty())
        if(not reaction_thermo_properties().helmholtz_energy.empty())
            return standardHelmholtzEnergyFromReaction(T, P, species, database, reaction_thermo_properties());

    if(hasThermoParamsHKF(species, database))
        return speciesThermoStateHKF(T, P, species, database).helmholtz_energy;

    Exception exception;
    exception.error << "Cannot calculate the standard Helmholtz energy of species " << species << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

auto standardInternalEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    const auto species_thermo_properties = getSpeciesThermoProperties(species, database);
    if(not species_thermo_properties.empty())
        if(not species_thermo_properties().internal_energy.empty())
            return ThermoScalar(species_thermo_properties().internal_energy(T, P), 0.0, 0.0);

    const auto reaction_thermo_properties = getReactionThermoProperties(species, database);
    if(not reaction_thermo_properties.empty())
        if(not reaction_thermo_properties().internal_energy.empty())
            return standardInternalEnergyFromReaction(T, P, species, database, reaction_thermo_properties());

    if(hasThermoParamsHKF(species, database))
        return speciesThermoStateHKF(T, P, species, database).internal_energy;

    Exception exception;
    exception.error << "Cannot calculate the standard internal energy of species " << species << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

auto standardEnthalpyEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    const auto species_thermo_properties = getSpeciesThermoProperties(species, database);
    if(not species_thermo_properties.empty())
        if(not species_thermo_properties().enthalpy.empty())
            return ThermoScalar(species_thermo_properties().enthalpy(T, P), 0.0, 0.0);

    const auto reaction_thermo_properties = getReactionThermoProperties(species, database);
    if(not reaction_thermo_properties.empty())
        if(not reaction_thermo_properties().enthalpy.empty())
            return standardEnthalpyFromReaction(T, P, species, database, reaction_thermo_properties());

    if(hasThermoParamsHKF(species, database))
        return speciesThermoStateHKF(T, P, species, database).enthalpy;

    Exception exception;
    exception.error << "Cannot calculate the standard enthalpy of species " << species << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

auto standardEntropy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    const auto species_thermo_properties = getSpeciesThermoProperties(species, database);
    if(not species_thermo_properties.empty())
        if(not species_thermo_properties().entropy.empty())
            return ThermoScalar(species_thermo_properties().entropy(T, P), 0.0, 0.0);

    const auto reaction_thermo_properties = getReactionThermoProperties(species, database);
    if(not reaction_thermo_properties.empty())
        if(not reaction_thermo_properties().entropy.empty())
            return standardEntropyFromReaction(T, P, species, database, reaction_thermo_properties());

    if(hasThermoParamsHKF(species, database))
        return speciesThermoStateHKF(T, P, species, database).entropy;

    Exception exception;
    exception.error << "Cannot calculate the standard entropy of species " << species << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

auto standardVolume(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    const auto species_thermo_properties = getSpeciesThermoProperties(species, database);
    if(not species_thermo_properties.empty())
        if(not species_thermo_properties().volume.empty())
            return ThermoScalar(species_thermo_properties().volume(T, P), 0.0, 0.0);

    const auto reaction_thermo_properties = getReactionThermoProperties(species, database);
    if(not reaction_thermo_properties.empty())
        if(not reaction_thermo_properties().volume.empty())
            return standardVolumeFromReaction(T, P, species, database, reaction_thermo_properties());

    if(hasThermoParamsHKF(species, database))
        return speciesThermoStateHKF(T, P, species, database).volume;

    Exception exception;
    exception.error << "Cannot calculate the standard volume of species " << species << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

auto standardHeatCapacity(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    const auto species_thermo_properties = getSpeciesThermoProperties(species, database);
    if(not species_thermo_properties.empty())
        if(not species_thermo_properties().heat_capacity.empty())
            return ThermoScalar(species_thermo_properties().heat_capacity(T, P), 0.0, 0.0);

    const auto reaction_thermo_properties = getReactionThermoProperties(species, database);
    if(not reaction_thermo_properties.empty())
        if(not reaction_thermo_properties().heat_capacity.empty())
            return standardHeatCapacityFromReaction(T, P, species, database, reaction_thermo_properties());

    if(hasThermoParamsHKF(species, database))
        return speciesThermoStateHKF(T, P, species, database).heat_capacity;

    Exception exception;
    exception.error << "Cannot calculate the standard heat capacity of species " << species << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

} // namespace Reaktor
