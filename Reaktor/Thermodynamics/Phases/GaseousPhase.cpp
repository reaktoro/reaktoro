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

#include "GaseousPhase.hpp"

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivityDuanSun.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivityIdeal.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivityPengRobinson.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivitySpycherPruess.hpp>
#include <Reaktor/Thermodynamics/Activity/GaseousActivitySpycherReed.hpp>

namespace Reaktor {

GaseousPhase::GaseousPhase()
: GaseousMixture()
{}

GaseousPhase::GaseousPhase(const std::vector<GaseousSpecies>& species)
: GaseousMixture(species), activities$(species.size())
{
    for(const auto& iter : species)
        setActivityModelIdeal(iter.name);
}

auto GaseousPhase::setActivityModel(const std::string& species, const GaseousActivity& activity) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        activities$[ispecies] = activity;
}

auto GaseousPhase::setActivityModelIdeal(const std::string& species) -> void
{
    const Index ispecies = indexSpecies(species);

    if(ispecies < numSpecies())
        activities$[ispecies] = gaseousActivityIdeal(species, *this);
}

auto GaseousPhase::setActivityModelDuanSunCO2() -> void
{
    const Index ispecies = indexSpecies("CO2(g)");

    if(ispecies < numSpecies())
        activities$[ispecies] = gaseousActivityDuanSunCO2(*this);
}

auto GaseousPhase::setActivityModelSpycherPruessH2OCO2() -> void
{
    const Index iH2O = indexSpecies("H2O(g)");
    const Index iCO2 = indexSpecies("CO2(g)");

    const auto& functions = gaseousActivitySpycherPruessH2OCO2(*this);

    if(iH2O < numSpecies()) activities$[iH2O] = functions[0];
    if(iCO2 < numSpecies()) activities$[iCO2] = functions[1];
}

auto GaseousPhase::setActivityModelSpycherReedH2OCO2CH4() -> void
{
    const Index iH2O = indexSpecies("H2O(g)");
    const Index iCO2 = indexSpecies("CO2(g)");
    const Index iCH4 = indexSpecies("CH4(g)");

    const auto& functions = gaseousActivitySpycherReedH2OCO2CH4(*this);

    if(iH2O < numSpecies()) activities$[iH2O] = functions[0];
    if(iCO2 < numSpecies()) activities$[iCO2] = functions[1];
    if(iCH4 < numSpecies()) activities$[iCH4] = functions[2];
}

auto GaseousPhase::setActivityModelPengRobinson(const std::string& species) -> void
{
    const Index idx_species = indexSpecies(species);
    if(idx_species < numSpecies())
        activities$[idx_species] = gaseousActivityPengRobinson(species, *this);
}

auto GaseousPhase::concentrations(const Vector& n) const -> Vector
{
    // The total amount of moles in the gaseous phase
    const double ntotal = n.sum();

    // Check if the phase has zero number of moles
    if(ntotal == 0.0) return zeros(n.rows());

    // Calculate the molar fractions of the gaseous species
    return n/ntotal;
}

auto GaseousPhase::activities(double T, double P, const Vector& n) const -> ChemicalVector
{
    GaseousMixtureState s = state(T, P, n);
    const unsigned N = numSpecies();
    ChemicalVector a(N, N);
    for(unsigned i = 0; i < N; ++i)
        a.row(i) = activities$[i](s);
    return a;
}

auto createPhase(const GaseousPhase& phase) -> Phase
{
//    // Create the gaseous species as Species instances
//    std::vector<Species> species;
//    for(const GaseousSpecies& iter : phase.species())
//        species.push_back(createSpecies(iter));
//
//    // Define the concentration function of the gaseous phase
//    Concentration concentration = [=](const Vector& n) -> Vector
//    {
//        return phase.concentrations(n);
//    };
//
//    // Define the activity function of the gaseous phase
//    Activity activity = [=](double T, double P, const Vector& n)
//    {
//        return phase.activities(T, P, n);
//    };
//
//    Phase converted;
//    converted.setName("Gaseous");
//    converted.setSpecies(species);
//    converted.setConcentration(concentration);
//    converted.setActivity(activity);
//
//    return converted;
}

} // namespace Reaktor
