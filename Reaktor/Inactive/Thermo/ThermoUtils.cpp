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

#include "ThermoUtils.hpp"

// C++ includes
#include <cmath>
#include <functional>
#include <sstream>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/Database.hpp>
#include <Reaktor/Core/ReactionEquation.hpp>
#include <Reaktor/External/Units/Units.hpp>
#include <Reaktor/Math/BilinearInterpolator.hpp>
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Species/GaseousSpecies.hpp>
#include <Reaktor/Species/MineralSpecies.hpp>
#include <Reaktor/Thermo/SpeciesThermo.hpp>
#include <Reaktor/Thermo/ThermoUtils.hpp>

namespace Reaktor {
namespace internal {

// The default temperature points for the interpolation of the thermodynamic quantities (in units of celsius)
const std::vector<double> tPoints = { 0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300 };

// The default pressure points for the interpolation of the thermodynamic quantities (in units of bar)
const std::vector<double> pPoints = { 1, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600 };

template<typename SpeciesType>
auto chemicalPotentialFn(const SpeciesType& species) -> BilinearInterpolator
{
    return chemicalPotentialFn(species, internal::tPoints, "degC", internal::pPoints, "bar");
}

template<typename SpeciesType>
auto chemicalPotentialFn(const SpeciesType& species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator
{
    // Convert the temperatures and pressures points to K and Pa respectively
    for(auto& x : tPoints) x = units::convert(x, tUnit, "K");
    for(auto& x : pPoints) x = units::convert(x, pUnit, "Pa");

    auto f = [&](double T, double P)
    {
        return speciesThermo(T, P, species).gibbs;
    };

    return BilinearInterpolator(tPoints, pPoints, f);
}

inline auto nonExistentSpeciesError(std::string species) -> void
{
    Exception exception;
    exception.error << "Unable to generate interpolation data for species " << species << ".";
    exception.reason << "This species is not contained in the database.";
    raise(exception);
}

} /* namespace internal */

auto chemicalPotentialFn(const Database& database, std::string species) -> BilinearInterpolator
{
    return chemicalPotentialFn(database, species, internal::tPoints, "degC", internal::pPoints, "bar");
}

auto chemicalPotentialFn(const Database& database, std::string species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator
{
    if(database.containsAqueousSpecies(species))
        return chemicalPotentialFn(database.aqueousSpecies(species), tPoints, tUnit, pPoints, pUnit);
    if(database.containsGaseousSpecies(species))
        return chemicalPotentialFn(database.gaseousSpecies(species), tPoints, tUnit, pPoints, pUnit);
    if(database.containsMineralSpecies(species))
        return chemicalPotentialFn(database.mineralSpecies(species), tPoints, tUnit, pPoints, pUnit);
    internal::nonExistentSpeciesError(species);
    return BilinearInterpolator();
}

auto chemicalPotentialFn(const AqueousSpecies& species) -> BilinearInterpolator
{
    return internal::chemicalPotentialFn(species);
}

auto chemicalPotentialFn(const AqueousSpecies& species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator
{
    return internal::chemicalPotentialFn(species, tPoints, tUnit, pPoints, pUnit);
}

auto chemicalPotentialFn(const GaseousSpecies& species) -> BilinearInterpolator
{
    return internal::chemicalPotentialFn(species);
}

auto chemicalPotentialFn(const GaseousSpecies& species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator
{
    return internal::chemicalPotentialFn(species, tPoints, tUnit, pPoints, pUnit);
}

auto chemicalPotentialFn(const MineralSpecies& species) -> BilinearInterpolator
{
    return internal::chemicalPotentialFn(species);
}

auto chemicalPotentialFn(const MineralSpecies& species,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator
{
    return internal::chemicalPotentialFn(species, tPoints, tUnit, pPoints, pUnit);
}

auto createEquilibriumConstant(const Database& database, const ReactionEquation& equation) -> BilinearInterpolator
{
    return createEquilibriumConstant(database, equation, internal::tPoints, "degC", internal::pPoints, "bar");
}

auto createEquilibriumConstant(const Database& database, const ReactionEquation& equation,
    std::vector<double> tPoints, std::string tUnit,
    std::vector<double> pPoints, std::string pUnit) -> BilinearInterpolator
{
    // The chemical potentials of the species in the reaction
    std::vector<BilinearInterpolator> potentials;

    // The stoichiometries of the species in the reaction
    std::vector<double> stoichiometries;

    for(const auto& pair : equation)
    {
        const auto& species = pair.first;
        const auto& stoichiometry = pair.second;

        potentials.push_back(chemicalPotentialFn(database, species, tPoints, tUnit, pPoints, pUnit));
        stoichiometries.push_back(stoichiometry);
    }

    auto f = [&](double T, double P)
    {
        double lnk = 0.0;
        for(unsigned i = 0; i < stoichiometries.size(); ++i)
            lnk += stoichiometries[i] * potentials[i](T, P);

        // The universal gas constant (in units of J/(mol*K))
        const double R = 8.3144621;

        return std::exp(-lnk/(R*T));
    };

    // Convert the temperatures and pressures points to K and Pa respectively
    for(auto& x : tPoints) x = units::convert(x, tUnit, "K");
    for(auto& x : pPoints) x = units::convert(x, pUnit, "Pa");

    return BilinearInterpolator(tPoints, pPoints, f);
}

//BilinearInterpolator
//chemicalPotentialFnFromLogK(const Database& database, const AqueousSpecies& species)
//{
//    const ThermoData::EquilibriumConstant& logk = species.thermoData().logk;
//
//    // Create the standard chemical potential function of each species
//    std::vector<ChemicalPotential> potentials;
//
//    for(const auto& pair : logk.reaction)
//        if(pair.first != species.name())
//        potentials.push_back(chemicalPotentialFn(database, pair.first));
//
//
//    ReactionEquation reaction = logk.reaction;
//
//    // Remove the entry in the reaction
//    std::remove_if(reaction.begin(), reaction.end(),
//        [&](const std::pair<std::string, double>& p)
//            { return p.first == species.name(); });
//}

} // namespace Reaktor
