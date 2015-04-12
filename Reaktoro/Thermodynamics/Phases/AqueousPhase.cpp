// Reaktoro is a C++ library for computational reaction modelling.
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

#include "AqueousPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityDrummond.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityDuanSun.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityHKF.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityIdeal.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityPitzer.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivityRumpf.hpp>
#include <Reaktoro/Thermodynamics/Activity/AqueousActivitySetschenow.hpp>
#include <Reaktoro/Thermodynamics/Water/WaterConstants.hpp>

namespace Reaktoro {

AqueousPhase::AqueousPhase()
: AqueousMixture()
{}

AqueousPhase::AqueousPhase(const std::vector<AqueousSpecies>& species)
: AqueousMixture(species), activity_fns(species.size())
{
    for(const auto& iter : species)
        setActivityModelSetschenow(iter.name(), 0.1);
}

auto AqueousPhase::setActivityModel(std::string species, const AqueousActivityFunction& activity) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        activity_fns[ispecies] = activity;
}

auto AqueousPhase::setActivityModelIdeal(std::string species) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        activity_fns[ispecies] = aqueousActivityIdeal(species, *this);
}

auto AqueousPhase::setActivityModelSetschenow(std::string species, double b) -> void
{
    const Index ispecies = indexSpecies(species);

    if(ispecies < numSpecies())
        activity_fns[ispecies] = aqueousActivitySetschenow(species, *this, b);
}

auto AqueousPhase::setActivityModelDuanSunCO2() -> void
{
    const Index ispecies = indexSpecies("CO2(aq)");

    if(ispecies < numSpecies())
        activity_fns[ispecies] = aqueousActivityDuanSunCO2(*this);
}

auto AqueousPhase::setActivityModelDrummondCO2() -> void
{
    const Index ispecies = indexSpecies("CO2(aq)");

    if(ispecies < numSpecies())
        activity_fns[ispecies] = aqueousActivityDrummondCO2(*this);
}

auto AqueousPhase::setActivityModelRumpfCO2() -> void
{
    const Index ispecies = indexSpecies("CO2(aq)");

    if(ispecies < numSpecies())
        activity_fns[ispecies] = aqueousActivityRumpfCO2(*this);
}

auto AqueousPhase::setActivityModelHKFWater() -> void
{
    const Index ispecies = indexSpecies("H2O(l)");

    if(ispecies < numSpecies())
        activity_fns[ispecies] = aqueousActivityHKFWater(*this);
}

auto AqueousPhase::setActivityModelHKFChargedSpecies() -> void
{
    for(Index idx : indicesChargedSpecies())
        activity_fns[idx] = aqueousActivityHKFCharged(species(idx).name(), *this);
}

auto AqueousPhase::setActivityModelPitzerWater() -> void
{
    const Index ispecies = indexSpecies("H2O(l)");

    if(ispecies < numSpecies())
        activity_fns[ispecies] = aqueousActivityPitzerWater(*this);
}

auto AqueousPhase::setActivityModelPitzerChargedSpecies() -> void
{
    for(Index idx : indicesChargedSpecies())
        activity_fns[idx] = aqueousActivityPitzerCharged(species(idx).name(), *this);
}

auto AqueousPhase::setActivityModelPitzerNeutralSpecies(std::string species) -> void
{
    const Index ispecies = indexSpecies(species);

    if(ispecies < numSpecies())
        activity_fns[ispecies] = aqueousActivityPitzerNeutral(species, *this);
}

auto AqueousPhase::concentrations(double T, double P, const Vector& n) const -> ChemicalVector
{
    // Calculate the molalities of the species
    ChemicalVector c = molalities(n);

    // Calculate the molar fractions of the species
    ChemicalVector x = molarFractions(n);

    // The index of the water species
    const Index iH2O = indexWater();

    // Set the concentration of water to its molar fraction
    c.row(iH2O) = x.row(iH2O);

    return c;
}

auto AqueousPhase::activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector
{
    return activities(T, P, n)/concentrations(T, P, n);
}

auto AqueousPhase::activities(double T, double P, const Vector& n) const -> ChemicalVector
{
    AqueousMixtureState mixture_state = state(T, P, n);
    const unsigned nspecies = numSpecies();
    ChemicalVector a(nspecies, nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        a.row(i) = activity_fns[i](mixture_state);
    return a;
}

} // namespace Reaktoro
