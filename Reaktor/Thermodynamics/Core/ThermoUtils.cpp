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
#include <Reaktor/Thermodynamics/Core/Database.hpp>
#include <Reaktor/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesThermoState.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesThermoStateHKF.hpp>

namespace Reaktor {
namespace internal {

auto errorNonExistentSpecies(const std::string& name) -> void
{
    Exception exception;
    exception.error << "Cannot get an instance of the species named " << name << " in the database.";
    exception.reason << "There is no such species in the database.";
    RaiseError(exception);
}

auto standardGibbsEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardHelmholtzEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardInternalEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardEnthalpy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardEntropy(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardVolume(double T, double P, std::string species, const Database& database) -> ThermoScalar;

auto standardHeatCapacityCp(double T, double P, std::string species, const Database& database) -> ThermoScalar;

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
    return standardPropertyFromReaction(T, P, species, database, reaction, internal::standardGibbsEnergy, eval);
}

auto standardHelmholtzEnergyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.helmholtz_energy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, internal::standardHelmholtzEnergy, eval);
}

auto standardInternalEnergyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.internal_energy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, internal::standardInternalEnergy, eval);
}

auto standardEnthalpyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.enthalpy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, internal::standardEnthalpy, eval);
}

auto standardEntropyFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.entropy(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, internal::standardEntropy, eval);
}

auto standardVolumeFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.volume(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, internal::standardVolume, eval);
}

auto standardHeatCapacityCpFromReaction(double T, double P, std::string species, const Database& database, const ReactionThermoProperties& reaction) -> ThermoScalar
{
    auto eval = [&]() { return reaction.heat_capacity_cp(T, P); };
    return standardPropertyFromReaction(T, P, species, database, reaction, internal::standardHeatCapacityCp, eval);
}

template<typename SpeciesType>
auto standardGibbsEnergyHelper(double T, double P, const SpeciesType& species, const Database& database) -> ThermoScalar
{
    if(not species.thermo.properties.empty())
        if(not species.thermo.properties().gibbs_energy.empty())
            return {species.thermo.properties().gibbs_energy(T, P), 0.0, 0.0};
    if(not species.thermo.reaction.empty())
        return standardGibbsEnergyFromReaction(T, P, species.name, database, species.thermo.reaction());
    if(not species.thermo.hkf.empty())
        return speciesThermoStateHKF(T, P, species).gibbs_energy;
    Exception exception;
    exception.error << "Cannot calculate the standard Gibbs energy of species " << species.name << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

template<typename SpeciesType>
auto standardHelmholtzEnergyHelper(double T, double P, const SpeciesType& species, const Database& database) -> ThermoScalar
{
    if(not species.thermo.properties.empty())
        if(not species.thermo.properties().helmholtz_energy.empty())
            return {species.thermo.properties().helmholtz_energy(T, P), 0.0, 0.0};
    if(not species.thermo.reaction.empty())
        return standardHelmholtzEnergyFromReaction(T, P, species.name, database, species.thermo.reaction());
    if(not species.thermo.hkf.empty())
        return speciesThermoStateHKF(T, P, species).helmholtz_energy;
    Exception exception;
    exception.error << "Cannot calculate the standard Helmholtz energy of species " << species.name << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

template<typename SpeciesType>
auto standardInternalEnergyHelper(double T, double P, const SpeciesType& species, const Database& database) -> ThermoScalar
{
    if(not species.thermo.properties.empty())
        if(not species.thermo.properties().internal_energy.empty())
            return {species.thermo.properties().internal_energy(T, P), 0.0, 0.0};
    if(not species.thermo.reaction.empty())
        return standardInternalEnergyFromReaction(T, P, species.name, database, species.thermo.reaction());
    if(not species.thermo.hkf.empty())
        return speciesThermoStateHKF(T, P, species).internal_energy;
    Exception exception;
    exception.error << "Cannot calculate the standard internal energy of species " << species.name << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

template<typename SpeciesType>
auto standardEnthalpyHelper(double T, double P, const SpeciesType& species, const Database& database) -> ThermoScalar
{
    if(not species.thermo.properties.empty())
        if(not species.thermo.properties().enthalpy.empty())
            return {species.thermo.properties().enthalpy(T, P), 0.0, 0.0};
    if(not species.thermo.reaction.empty())
        return standardEnthalpyFromReaction(T, P, species.name, database, species.thermo.reaction());
    if(not species.thermo.hkf.empty())
        return speciesThermoStateHKF(T, P, species).enthalpy;
    Exception exception;
    exception.error << "Cannot calculate the standard enthalpy of species " << species.name << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

template<typename SpeciesType>
auto standardEntropyHelper(double T, double P, const SpeciesType& species, const Database& database) -> ThermoScalar
{
    if(not species.thermo.properties.empty())
        if(not species.thermo.properties().entropy.empty())
            return {species.thermo.properties().entropy(T, P), 0.0, 0.0};
    if(not species.thermo.reaction.empty())
        return standardEntropyFromReaction(T, P, species.name, database, species.thermo.reaction());
    if(not species.thermo.hkf.empty())
        return speciesThermoStateHKF(T, P, species).entropy;
    Exception exception;
    exception.error << "Cannot calculate the standard entropy of species " << species.name << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

template<typename SpeciesType>
auto standardVolumeHelper(double T, double P, const SpeciesType& species, const Database& database) -> ThermoScalar
{
    if(not species.thermo.properties.empty())
        if(not species.thermo.properties().volume.empty())
            return {species.thermo.properties().volume(T, P), 0.0, 0.0};
    if(not species.thermo.reaction.empty())
        return standardVolumeFromReaction(T, P, species.name, database, species.thermo.reaction());
    if(not species.thermo.hkf.empty())
        return speciesThermoStateHKF(T, P, species).volume;
    Exception exception;
    exception.error << "Cannot calculate the standard volume of species " << species.name << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

template<typename SpeciesType>
auto standardHeatCapacityCpHelper(double T, double P, const SpeciesType& species, const Database& database) -> ThermoScalar
{
    if(not species.thermo.properties.empty())
        if(not species.thermo.properties().heat_capacity_cp.empty())
            return {species.thermo.properties().heat_capacity_cp(T, P), 0.0, 0.0};
    if(not species.thermo.reaction.empty())
        return standardHeatCapacityCpFromReaction(T, P, species.name, database, species.thermo.reaction());
    if(not species.thermo.hkf.empty())
        return speciesThermoStateHKF(T, P, species).heat_capacity_cp;
    Exception exception;
    exception.error << "Cannot calculate the standard heat capacity of species " << species.name << ".";
    exception.reason << "The species instance has no thermodynamic data for such calculation.";
    RaiseError(exception);
    return {};
}

auto standardGibbsEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    if(database.containsAqueousSpecies(species))
        return standardGibbsEnergyHelper(T, P, database.aqueousSpecies(species), database);
    if(database.containsGaseousSpecies(species))
        return standardGibbsEnergyHelper(T, P, database.gaseousSpecies(species), database);
    if(database.containsMineralSpecies(species))
        return standardGibbsEnergyHelper(T, P, database.mineralSpecies(species), database);
    internal::errorNonExistentSpecies(species);
    return {};
}

auto standardHelmholtzEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    if(database.containsAqueousSpecies(species))
        return standardHelmholtzEnergyHelper(T, P, database.aqueousSpecies(species), database);
    if(database.containsGaseousSpecies(species))
        return standardHelmholtzEnergyHelper(T, P, database.gaseousSpecies(species), database);
    if(database.containsMineralSpecies(species))
        return standardHelmholtzEnergyHelper(T, P, database.mineralSpecies(species), database);
    internal::errorNonExistentSpecies(species);
    return {};
}

auto standardInternalEnergy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    if(database.containsAqueousSpecies(species))
        return standardInternalEnergyHelper(T, P, database.aqueousSpecies(species), database);
    if(database.containsGaseousSpecies(species))
        return standardInternalEnergyHelper(T, P, database.gaseousSpecies(species), database);
    if(database.containsMineralSpecies(species))
        return standardInternalEnergyHelper(T, P, database.mineralSpecies(species), database);
    internal::errorNonExistentSpecies(species);
    return {};
}

auto standardEnthalpy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    if(database.containsAqueousSpecies(species))
        return standardEnthalpyHelper(T, P, database.aqueousSpecies(species), database);
    if(database.containsGaseousSpecies(species))
        return standardEnthalpyHelper(T, P, database.gaseousSpecies(species), database);
    if(database.containsMineralSpecies(species))
        return standardEnthalpyHelper(T, P, database.mineralSpecies(species), database);
    internal::errorNonExistentSpecies(species);
    return {};
}

auto standardEntropy(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    if(database.containsAqueousSpecies(species))
        return standardEntropyHelper(T, P, database.aqueousSpecies(species), database);
    if(database.containsGaseousSpecies(species))
        return standardEntropyHelper(T, P, database.gaseousSpecies(species), database);
    if(database.containsMineralSpecies(species))
        return standardEntropyHelper(T, P, database.mineralSpecies(species), database);
    internal::errorNonExistentSpecies(species);
    return {};
}

auto standardVolume(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    if(database.containsAqueousSpecies(species))
        return standardVolumeHelper(T, P, database.aqueousSpecies(species), database);
    if(database.containsGaseousSpecies(species))
        return standardVolumeHelper(T, P, database.gaseousSpecies(species), database);
    if(database.containsMineralSpecies(species))
        return standardVolumeHelper(T, P, database.mineralSpecies(species), database);
    internal::errorNonExistentSpecies(species);
    return {};
}

auto standardHeatCapacityCp(double T, double P, std::string species, const Database& database) -> ThermoScalar
{
    if(database.containsAqueousSpecies(species))
        return standardHeatCapacityCpHelper(T, P, database.aqueousSpecies(species), database);
    if(database.containsGaseousSpecies(species))
        return standardHeatCapacityCpHelper(T, P, database.gaseousSpecies(species), database);
    if(database.containsMineralSpecies(species))
        return standardHeatCapacityCpHelper(T, P, database.mineralSpecies(species), database);
    internal::errorNonExistentSpecies(species);
    return {};
}









template<typename SpeciesType, typename PropertyFunction>
auto interpolatePropertyFromReaction(
    const std::vector<double>& temperatures,
    const std::vector<double>& pressures,
    const SpeciesType& species,
    const Database& database, PropertyFunction property) -> BilinearInterpolator
{
    auto func = [&](double T, double P)
    {
        return property(T, P, species.name, database, species.thermo.reaction());
    };
    return BilinearInterpolator(temperatures, pressures, func);
}

template<typename SpeciesType>
auto interpolateStandardGibbsEnergyFromReaction(const SpeciesType& species, const Database& database) -> BilinearInterpolator
{
    const auto& temperatures = species.thermo.reaction().gibbs_energy.xCoodinates();
    const auto& pressures = species.thermo.reaction().gibbs_energy.yCoodinates();
    return interpolatePropertyFromReaction(temperatures, pressures, species, database, internal::standardGibbsEnergyFromReaction);
}

template<typename SpeciesType>
auto interpolateStandardHelmholtzEnergyFromReaction(const SpeciesType& species, const Database& database) -> BilinearInterpolator
{
    const auto& temperatures = species.thermo.reaction().gibbs_energy.xCoodinates();
    const auto& pressures = species.thermo.reaction().gibbs_energy.yCoodinates();
    return interpolatePropertyFromReaction(temperatures, pressures, species, database, internal::standardHelmholtzEnergyFromReaction);
}

template<typename SpeciesType>
auto interpolateStandardInternalEnergyFromReaction(const SpeciesType& species, const Database& database) -> BilinearInterpolator
{
    const auto& temperatures = species.thermo.reaction().gibbs_energy.xCoodinates();
    const auto& pressures = species.thermo.reaction().gibbs_energy.yCoodinates();
    return interpolatePropertyFromReaction(temperatures, pressures, species, database, internal::standardInternalEnergyFromReaction);
}

template<typename SpeciesType>
auto interpolateStandardEnthalpyFromReaction(const SpeciesType& species, const Database& database) -> BilinearInterpolator
{
    const auto& temperatures = species.thermo.reaction().gibbs_energy.xCoodinates();
    const auto& pressures = species.thermo.reaction().gibbs_energy.yCoodinates();
    return interpolatePropertyFromReaction(temperatures, pressures, species, database, internal::standardEnthalpyFromReaction);
}

template<typename SpeciesType>
auto interpolateStandardEntropyFromReaction(const SpeciesType& species, const Database& database) -> BilinearInterpolator
{
    const auto& temperatures = species.thermo.reaction().gibbs_energy.xCoodinates();
    const auto& pressures = species.thermo.reaction().gibbs_energy.yCoodinates();
    return interpolatePropertyFromReaction(temperatures, pressures, species, database, internal::standardEntropyFromReaction);
}

template<typename SpeciesType>
auto interpolateStandardVolumeFromReaction(const SpeciesType& species, const Database& database) -> BilinearInterpolator
{
    const auto& temperatures = species.thermo.reaction().gibbs_energy.xCoodinates();
    const auto& pressures = species.thermo.reaction().gibbs_energy.yCoodinates();
    return interpolatePropertyFromReaction(temperatures, pressures, species, database, internal::standardVolumeFromReaction);
}

} // namespace internal

auto standardGibbsEnergy(double T, double P, const AqueousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardGibbsEnergyHelper(T, P, species, database);
}

auto standardGibbsEnergy(double T, double P, const GaseousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardGibbsEnergyHelper(T, P, species, database);
}

auto standardGibbsEnergy(double T, double P, const MineralSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardGibbsEnergyHelper(T, P, species, database);
}

auto standardHelmholtzEnergy(double T, double P, const AqueousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardHelmholtzEnergyHelper(T, P, species, database);
}

auto standardHelmholtzEnergy(double T, double P, const GaseousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardHelmholtzEnergyHelper(T, P, species, database);
}

auto standardHelmholtzEnergy(double T, double P, const MineralSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardHelmholtzEnergyHelper(T, P, species, database);
}

auto standardInternalEnergy(double T, double P, const AqueousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardInternalEnergyHelper(T, P, species, database);
}

auto standardInternalEnergy(double T, double P, const GaseousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardInternalEnergyHelper(T, P, species, database);
}

auto standardInternalEnergy(double T, double P, const MineralSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardInternalEnergyHelper(T, P, species, database);
}

auto standardEnthalpy(double T, double P, const AqueousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardEnthalpyHelper(T, P, species, database);
}

auto standardEnthalpy(double T, double P, const GaseousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardEnthalpyHelper(T, P, species, database);
}

auto standardEnthalpy(double T, double P, const MineralSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardEnthalpyHelper(T, P, species, database);
}

auto standardEntropy(double T, double P, const AqueousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardEntropyHelper(T, P, species, database);
}

auto standardEntropy(double T, double P, const GaseousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardEntropyHelper(T, P, species, database);
}

auto standardEntropy(double T, double P, const MineralSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardEntropyHelper(T, P, species, database);
}

auto standardVolume(double T, double P, const AqueousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardVolumeHelper(T, P, species, database);
}

auto standardVolume(double T, double P, const GaseousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardVolumeHelper(T, P, species, database);
}

auto standardVolume(double T, double P, const MineralSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardVolumeHelper(T, P, species, database);
}

auto standardHeatCapacityCp(double T, double P, const AqueousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardHeatCapacityCpHelper(T, P, species, database);
}

auto standardHeatCapacityCp(double T, double P, const GaseousSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardHeatCapacityCpHelper(T, P, species, database);
}

auto standardHeatCapacityCp(double T, double P, const MineralSpecies& species, const Database& database) -> ThermoScalar
{
    return internal::standardHeatCapacityCpHelper(T, P, species, database);
}

} // namespace Reaktor
