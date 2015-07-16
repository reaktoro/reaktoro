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

#include "GaseousPhase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Activity/GaseousActivityDuanSun.hpp>
#include <Reaktoro/Thermodynamics/Activity/GaseousActivityIdeal.hpp>
#include <Reaktoro/Thermodynamics/Activity/GaseousActivityPengRobinson.hpp>
#include <Reaktoro/Thermodynamics/Activity/GaseousActivitySpycherPruess.hpp>
#include <Reaktoro/Thermodynamics/Activity/GaseousActivitySpycherReed.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/GaseousChemicalModelIdeal.hpp>
#include <Reaktoro/Thermodynamics/Models/GaseousChemicalModelSpycherEtAl2003.hpp>

namespace Reaktoro {

struct GaseousPhase::Impl
{
    /// The gaseous mixture instance
    GaseousMixture mixture;

    /// The functions that calculates the activities of selected species
    std::map<Index, GaseousActivityFunction> activities;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance
    Impl(const GaseousMixture& mixture)
    : mixture(mixture)
    {}
};

GaseousPhase::GaseousPhase()
: Phase()
{}

GaseousPhase::GaseousPhase(const GaseousPhase& other)
: Phase(other), pimpl(new Impl(*other.pimpl))
{}

GaseousPhase::GaseousPhase(const GaseousMixture& mixture)
: pimpl(new Impl(mixture))
{
    // Convert the GaseousSpecies instances to Species instances
    std::vector<Species> species;
    for(const GaseousSpecies& x : mixture.species())
        species.push_back(x);

    // Set the Phase attributes
    setName("Gaseous");
    setSpecies(species);
    setReferenceState(PhaseReferenceState::IdealGas);
    setChemicalModelIdeal();
}

GaseousPhase::~GaseousPhase()
{}

auto GaseousPhase::operator=(GaseousPhase other) -> GaseousPhase&
{
    Phase::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto GaseousPhase::setChemicalModelIdeal() -> void
{
    // Create the gaseous chemical model
    PhaseChemicalModel model = gaseousChemicalModelIdeal(mixture());

    setChemicalModel(model);
}

auto GaseousPhase::setChemicalModelSpycherEtAl2003() -> void
{
    // Create the gaseous chemical model
    PhaseChemicalModel model = gaseousChemicalModelSpycherEtAl2003(mixture());

    setChemicalModel(model);
}

auto GaseousPhase::setActivityModel(std::string species, const GaseousActivityFunction& activity) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = activity;
}

auto GaseousPhase::setActivityModelIdeal(std::string species) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = gaseousActivityIdeal(species, mixture());
}

auto GaseousPhase::setActivityModelDuanSunCO2() -> void
{
    const Index ispecies = indexSpecies("CO2(g)");
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = gaseousActivityDuanSunCO2(mixture());
}

auto GaseousPhase::setActivityModelSpycherPruessH2OCO2() -> void
{
    const Index iH2O = indexSpecies("H2O(g)");
    const Index iCO2 = indexSpecies("CO2(g)");

    const auto& functions = gaseousActivitySpycherPruessH2OCO2(mixture());

    if(iH2O < numSpecies()) pimpl->activities[iH2O] = functions[0];
    if(iCO2 < numSpecies()) pimpl->activities[iCO2] = functions[1];
}

auto GaseousPhase::setActivityModelSpycherReedH2OCO2CH4() -> void
{
    const Index iH2O = indexSpecies("H2O(g)");
    const Index iCO2 = indexSpecies("CO2(g)");
    const Index iCH4 = indexSpecies("CH4(g)");

    const auto& functions = gaseousActivitySpycherReedH2OCO2CH4(mixture());

    if(iH2O < numSpecies()) pimpl->activities[iH2O] = functions[0];
    if(iCO2 < numSpecies()) pimpl->activities[iCO2] = functions[1];
    if(iCH4 < numSpecies()) pimpl->activities[iCH4] = functions[2];
}

auto GaseousPhase::setActivityModelPengRobinson(std::string species) -> void
{
    const Index idx_species = indexSpecies(species);
    if(idx_species < numSpecies())
        pimpl->activities[idx_species] = gaseousActivityPengRobinson(species, mixture());
}

auto GaseousPhase::mixture() const -> const GaseousMixture&
{
    return pimpl->mixture;
}

} // namespace Reaktoro
