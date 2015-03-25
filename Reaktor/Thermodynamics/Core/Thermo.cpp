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

#include "Thermo.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/OptimizationUtils.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Thermodynamics/Core/Database.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesElectroState.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesElectroStateHKF.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesThermoState.hpp>
#include <Reaktor/Thermodynamics/Models/SpeciesThermoStateHKF.hpp>
#include <Reaktor/Thermodynamics/Species/AqueousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/GaseousSpecies.hpp>
#include <Reaktor/Thermodynamics/Species/MineralSpecies.hpp>
#include <Reaktor/Thermodynamics/Water/WaterElectroState.hpp>
#include <Reaktor/Thermodynamics/Water/WaterElectroStateJohnsonNorton.hpp>
#include <Reaktor/Thermodynamics/Water/WaterThermoState.hpp>
#include <Reaktor/Thermodynamics/Water/WaterThermoStateUtils.hpp>

namespace Reaktor {
namespace {

/// The signature of a function that calculates the thermodynamic state of water
using WaterThermoStateFunction =
    std::function<WaterThermoState(double, double)>;

/// The signature of a function that calculates the thermodynamic state of a species
using SpeciesThermoStateFunction =
    std::function<SpeciesThermoState(double, double, std::string)>;

/// The signature of a function that calculates the electrostatic state of water
using WaterElectroStateFunction =
    std::function<WaterElectroState(double, double)>;

auto errorNonExistentSpecies(const std::string& name) -> void
{
    Exception exception;
    exception.error << "Cannot get an instance of the species `" << name << "` in the database.";
    exception.reason << "There is no such species in the database.";
    RaiseError(exception);
}

} // namespace

struct Thermo::Impl
{
    /// The database instance
    Database database;

    /// The Haar-Gallagher-Kell (1984) equation of state for water
    WaterThermoStateFunction water_thermo_state_hgk_fn;

    /// The Wagner and Pruss (1995) equation of state for water
    WaterThermoStateFunction water_thermo_state_wagner_pruss_fn;

    /// The Johnson and Norton equation of state for the electrostatic state of water
    WaterElectroStateFunction water_eletro_state_fn;

    /// The HKF equation of state for the thermodynamic state of aqueous, gaseous and mineral species
    SpeciesThermoStateFunction species_thermo_state_hkf_fn;

    Impl()
    {}

    Impl(const Database& database)
    : database(database)
    {
        // Initialize the Haar-Gallagher-Kell (1984) equation of state for water
        water_thermo_state_hgk_fn = [](double T, double P)
        {
            return Reaktor::waterThermoStateHGK(T, P);
        };

        water_thermo_state_hgk_fn = memoize(water_thermo_state_hgk_fn);

        // Initialize the Wagner and Pruss (1995) equation of state for water
        water_thermo_state_wagner_pruss_fn = [](double T, double P)
        {
            return Reaktor::waterThermoStateWagnerPruss(T, P);
        };

        water_thermo_state_wagner_pruss_fn = memoize(water_thermo_state_wagner_pruss_fn);

        // Initialize the Johnson and Norton equation of state for the electrostatic state of water
        water_eletro_state_fn = [=](double T, double P)
        {
            const WaterThermoState wts = water_thermo_state_wagner_pruss_fn(T, P);
            return waterElectroStateJohnsonNorton(T, P, wts);
        };

        water_eletro_state_fn = memoize(water_eletro_state_fn);

        // Initialize the HKF equation of state for the thermodynamic state of aqueous, gaseous and mineral species
        species_thermo_state_hkf_fn = [=](double T, double P, std::string species)
        {
            return speciesThermoStateHKF(T, P, species);
        };

        species_thermo_state_hkf_fn = memoize(species_thermo_state_hkf_fn);
    }

    auto speciesThermoStateHKF(double T, double P, std::string species) -> SpeciesThermoState
    {
        if(database.containsAqueousSpecies(species))
            return aqueousSpeciesThermoStateHKF(T, P, database.aqueousSpecies(species));
        if(database.containsGaseousSpecies(species))
            return Reaktor::speciesThermoStateHKF(T, P, database.gaseousSpecies(species));
        if(database.containsMineralSpecies(species))
            return Reaktor::speciesThermoStateHKF(T, P, database.mineralSpecies(species));
        errorNonExistentSpecies(species);
        return {};
    }

    auto aqueousSpeciesThermoStateHKF(double T, double P, const AqueousSpecies& species) -> SpeciesThermoState
    {
        const WaterThermoState wts = water_thermo_state_wagner_pruss_fn(T, P);

        if(species.name() == "H2O(l)")
            return speciesThermoStateSolventHKF(T, P, wts);

        const WaterElectroState wes = water_eletro_state_fn(T, P);

        FunctionG g = functionG(T, P, wts);

        SpeciesElectroState aes = speciesElectroStateHKF(g, species);

        return speciesThermoStateSoluteHKF(T, P, species, aes, wes);
    }

    auto standardGibbsEnergy(double T, double P, std::string species) -> ThermoScalar
    {
        const auto species_thermo_properties = getSpeciesThermoProperties(species);
        if(not species_thermo_properties.empty())
            if(not species_thermo_properties().gibbs_energy.empty())
                return ThermoScalar(species_thermo_properties().gibbs_energy(T, P), 0.0, 0.0);

        const auto reaction_thermo_properties = getReactionThermoProperties(species);
        if(not reaction_thermo_properties.empty())
            if(not reaction_thermo_properties.get().gibbs_energy.empty())
                return standardGibbsEnergyFromReaction(T, P, species, reaction_thermo_properties.get());

        if(hasThermoParamsHKF(species))
            return species_thermo_state_hkf_fn(T, P, species).gibbs_energy;

        Exception exception;
        exception.error << "Cannot calculate the standard Gibbs energy of species `" << species << "`.";
        exception.reason << "The species instance has no thermodynamic data for such calculation.";
        RaiseError(exception);
        return {};
    }

    auto standardHelmholtzEnergy(double T, double P, std::string species) -> ThermoScalar
    {
        const auto species_thermo_properties = getSpeciesThermoProperties(species);
        if(not species_thermo_properties.empty())
            if(not species_thermo_properties().helmholtz_energy.empty())
                return ThermoScalar(species_thermo_properties().helmholtz_energy(T, P), 0.0, 0.0);

        const auto reaction_thermo_properties = getReactionThermoProperties(species);
        if(not reaction_thermo_properties.empty())
            if(not reaction_thermo_properties.get().helmholtz_energy.empty())
                return standardHelmholtzEnergyFromReaction(T, P, species, reaction_thermo_properties.get());

        if(hasThermoParamsHKF(species))
            return species_thermo_state_hkf_fn(T, P, species).helmholtz_energy;

        Exception exception;
        exception.error << "Cannot calculate the standard Helmholtz energy of species `" << species << "`.";
        exception.reason << "The species instance has no thermodynamic data for such calculation.";
        RaiseError(exception);
        return {};
    }

    auto standardInternalEnergy(double T, double P, std::string species) -> ThermoScalar
    {
        const auto species_thermo_properties = getSpeciesThermoProperties(species);
        if(not species_thermo_properties.empty())
            if(not species_thermo_properties().internal_energy.empty())
                return ThermoScalar(species_thermo_properties().internal_energy(T, P), 0.0, 0.0);

        const auto reaction_thermo_properties = getReactionThermoProperties(species);
        if(not reaction_thermo_properties.empty())
            if(not reaction_thermo_properties.get().internal_energy.empty())
                return standardInternalEnergyFromReaction(T, P, species, reaction_thermo_properties.get());

        if(hasThermoParamsHKF(species))
            return species_thermo_state_hkf_fn(T, P, species).internal_energy;

        Exception exception;
        exception.error << "Cannot calculate the standard internal energy of species `" << species << "`.";
        exception.reason << "The species instance has no thermodynamic data for such calculation.";
        RaiseError(exception);
        return {};
    }

    auto standardEnthalpy(double T, double P, std::string species) -> ThermoScalar
    {
        const auto species_thermo_properties = getSpeciesThermoProperties(species);
        if(not species_thermo_properties.empty())
            if(not species_thermo_properties().enthalpy.empty())
                return ThermoScalar(species_thermo_properties().enthalpy(T, P), 0.0, 0.0);

        const auto reaction_thermo_properties = getReactionThermoProperties(species);
        if(not reaction_thermo_properties.empty())
            if(not reaction_thermo_properties.get().enthalpy.empty())
                return standardEnthalpyFromReaction(T, P, species, reaction_thermo_properties.get());

        if(hasThermoParamsHKF(species))
            return species_thermo_state_hkf_fn(T, P, species).enthalpy;

        Exception exception;
        exception.error << "Cannot calculate the standard enthalpy of species `" << species << "`.";
        exception.reason << "The species instance has no thermodynamic data for such calculation.";
        RaiseError(exception);
        return {};
    }

    auto standardEntropy(double T, double P, std::string species) -> ThermoScalar
    {
        const auto species_thermo_properties = getSpeciesThermoProperties(species);
        if(not species_thermo_properties.empty())
            if(not species_thermo_properties().entropy.empty())
                return ThermoScalar(species_thermo_properties().entropy(T, P), 0.0, 0.0);

        const auto reaction_thermo_properties = getReactionThermoProperties(species);
        if(not reaction_thermo_properties.empty())
            if(not reaction_thermo_properties.get().entropy.empty())
                return standardEntropyFromReaction(T, P, species, reaction_thermo_properties.get());

        if(hasThermoParamsHKF(species))
            return species_thermo_state_hkf_fn(T, P, species).entropy;

        Exception exception;
        exception.error << "Cannot calculate the standard entropy of species `" << species << "`.";
        exception.reason << "The species instance has no thermodynamic data for such calculation.";
        RaiseError(exception);
        return {};
    }

    auto standardVolume(double T, double P, std::string species) -> ThermoScalar
    {
        const auto species_thermo_properties = getSpeciesThermoProperties(species);
        if(not species_thermo_properties.empty())
            if(not species_thermo_properties().volume.empty())
                return ThermoScalar(species_thermo_properties().volume(T, P), 0.0, 0.0);

        const auto reaction_thermo_properties = getReactionThermoProperties(species);
        if(not reaction_thermo_properties.empty())
            if(not reaction_thermo_properties.get().volume.empty())
                return standardVolumeFromReaction(T, P, species, reaction_thermo_properties.get());

        if(hasThermoParamsHKF(species))
            return species_thermo_state_hkf_fn(T, P, species).volume;

        Exception exception;
        exception.error << "Cannot calculate the standard volume of species `" << species << "`.";
        exception.reason << "The species instance has no thermodynamic data for such calculation.";
        RaiseError(exception);
        return {};
    }

    auto standardHeatCapacity(double T, double P, std::string species) -> ThermoScalar
    {
        const auto species_thermo_properties = getSpeciesThermoProperties(species);
        if(not species_thermo_properties.empty())
            if(not species_thermo_properties().heat_capacity.empty())
                return ThermoScalar(species_thermo_properties().heat_capacity(T, P), 0.0, 0.0);

        const auto reaction_thermo_properties = getReactionThermoProperties(species);
        if(not reaction_thermo_properties.empty())
            if(not reaction_thermo_properties.get().heat_capacity.empty())
                return standardHeatCapacityFromReaction(T, P, species, reaction_thermo_properties.get());

        if(hasThermoParamsHKF(species))
            return species_thermo_state_hkf_fn(T, P, species).heat_capacity;

        Exception exception;
        exception.error << "Cannot calculate the standard heat capacity of species `" << species << "`.";
        exception.reason << "The species instance has no thermodynamic data for such calculation.";
        RaiseError(exception);
        return {};
    }

    auto getSpeciesThermoProperties(std::string species) -> Optional<SpeciesThermoProperties>
    {
        if(database.containsAqueousSpecies(species))
            return database.aqueousSpecies(species).thermoData().properties;
        if(database.containsGaseousSpecies(species))
            return database.gaseousSpecies(species).thermoData().properties;
        if(database.containsMineralSpecies(species))
            return database.mineralSpecies(species).thermoData().properties;
        errorNonExistentSpecies(species);
        return {};
    }

    auto getReactionThermoProperties(std::string species) -> Optional<ReactionThermoProperties>
    {
        if(database.containsAqueousSpecies(species))
            return database.aqueousSpecies(species).thermoData().reaction;
        if(database.containsGaseousSpecies(species))
            return database.gaseousSpecies(species).thermoData().reaction;
        if(database.containsMineralSpecies(species))
            return database.mineralSpecies(species).thermoData().reaction;
        errorNonExistentSpecies(species);
        return {};
    }

    auto hasThermoParamsHKF(std::string species) -> bool
    {
        if(species == "H2O(l)") return true;
        if(database.containsAqueousSpecies(species))
            return not database.aqueousSpecies(species).thermoData().hkf.empty();
        if(database.containsGaseousSpecies(species))
            return not database.gaseousSpecies(species).thermoData().hkf.empty();
        if(database.containsMineralSpecies(species))
            return not database.mineralSpecies(species).thermoData().hkf.empty();
        errorNonExistentSpecies(species);
        return {};
    }

    template<typename PropertyFunction, typename EvalFunction>
    auto standardPropertyFromReaction(double T, double P, std::string species,
        const ReactionThermoProperties& reaction, PropertyFunction property,
        EvalFunction eval) -> ThermoScalar
    {
        const double stoichiometry = reaction.equation.stoichiometry(species);

        Assert(stoichiometry, "Cannot calculate the thermodynamic property of "
            "species `" + species + "` using its reaction data.", "This species "
            "is not present in the reaction equation `" +
            std::string(reaction.equation) + "` or has zero stoichiometry.");

        double sum = 0.0;
        for(auto pair : reaction.equation)
        {
            const auto reactant = pair.first;
            const auto stoichiometry = pair.second;
            if(reactant != species)
                sum -= stoichiometry * property(T, P, reactant).val;
        }
        sum += eval();
        sum /= stoichiometry;

        return {sum, 0.0, 0.0};
    }

    auto standardGibbsEnergyFromReaction(double T, double P, std::string species, const ReactionThermoProperties& reaction) -> ThermoScalar
    {
        auto eval = [&]() { return reaction.gibbs_energy(T, P); };
        auto property = std::bind(&Impl::standardGibbsEnergy, *this, _1, _2, _3);
        return standardPropertyFromReaction(T, P, species, reaction, property, eval);
    }

    auto standardHelmholtzEnergyFromReaction(double T, double P, std::string species, const ReactionThermoProperties& reaction) -> ThermoScalar
    {
        auto eval = [&]() { return reaction.helmholtz_energy(T, P); };
        auto property = std::bind(&Impl::standardHelmholtzEnergy, *this, _1, _2, _3);
        return standardPropertyFromReaction(T, P, species, reaction, property, eval);
    }

    auto standardInternalEnergyFromReaction(double T, double P, std::string species, const ReactionThermoProperties& reaction) -> ThermoScalar
    {
        auto eval = [&]() { return reaction.internal_energy(T, P); };
        auto property = std::bind(&Impl::standardInternalEnergy, *this, _1, _2, _3);
        return standardPropertyFromReaction(T, P, species, reaction, property, eval);
    }

    auto standardEnthalpyFromReaction(double T, double P, std::string species, const ReactionThermoProperties& reaction) -> ThermoScalar
    {
        auto eval = [&]() { return reaction.enthalpy(T, P); };
        auto property = std::bind(&Impl::standardEnthalpy, *this, _1, _2, _3);
        return standardPropertyFromReaction(T, P, species, reaction, property, eval);
    }

    auto standardEntropyFromReaction(double T, double P, std::string species, const ReactionThermoProperties& reaction) -> ThermoScalar
    {
        auto eval = [&]() { return reaction.entropy(T, P); };
        auto property = std::bind(&Impl::standardEntropy, *this, _1, _2, _3);
        return standardPropertyFromReaction(T, P, species, reaction, property, eval);
    }

    auto standardVolumeFromReaction(double T, double P, std::string species, const ReactionThermoProperties& reaction) -> ThermoScalar
    {
        auto eval = [&]() { return reaction.volume(T, P); };
        auto property = std::bind(&Impl::standardVolume, *this, _1, _2, _3);
        return standardPropertyFromReaction(T, P, species, reaction, property, eval);
    }

    auto standardHeatCapacityFromReaction(double T, double P, std::string species, const ReactionThermoProperties& reaction) -> ThermoScalar
    {
        auto eval = [&]() { return reaction.heat_capacity(T, P); };
        auto property = std::bind(&Impl::standardHeatCapacity, *this, _1, _2, _3);
        return standardPropertyFromReaction(T, P, species, reaction, property, eval);
    }
};

Thermo::Thermo(const Database& database)
: pimpl(new Impl(database))
{}

auto Thermo::standardGibbsEnergy(double T, double P, std::string species) -> ThermoScalar
{
    return pimpl->standardGibbsEnergy(T, P, species);
}

auto Thermo::standardHelmholtzEnergy(double T, double P, std::string species) -> ThermoScalar
{
    return pimpl->standardHelmholtzEnergy(T, P, species);
}

auto Thermo::standardInternalEnergy(double T, double P, std::string species) -> ThermoScalar
{
    return pimpl->standardInternalEnergy(T, P, species);
}

auto Thermo::standardEnthalpy(double T, double P, std::string species) -> ThermoScalar
{
    return pimpl->standardEnthalpy(T, P, species);
}

auto Thermo::standardEntropy(double T, double P, std::string species) -> ThermoScalar
{
    return pimpl->standardEntropy(T, P, species);
}

auto Thermo::standardVolume(double T, double P, std::string species) -> ThermoScalar
{
    return pimpl->standardVolume(T, P, species);
}

auto Thermo::standardHeatCapacity(double T, double P, std::string species) -> ThermoScalar
{
    return pimpl->standardHeatCapacity(T, P, species);
}

auto Thermo::speciesThermoStateHKF(double T, double P, std::string species) -> SpeciesThermoState
{
    return pimpl->species_thermo_state_hkf_fn(T, P, species);
}

auto Thermo::waterThermoStateHGK(double T, double P) -> WaterThermoState
{
    return pimpl->water_thermo_state_hgk_fn(T, P);
}

auto Thermo::waterThermoStateWagnerPruss(double T, double P) -> WaterThermoState
{
    return pimpl->water_thermo_state_wagner_pruss_fn(T, P);
}

} // namespace Reaktor
