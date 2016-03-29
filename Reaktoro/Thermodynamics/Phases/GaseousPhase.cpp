// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <Reaktoro/Thermodynamics/Mixtures/GaseousMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/GaseousChemicalModelCubicEOS.hpp>
#include <Reaktoro/Thermodynamics/Models/GaseousChemicalModelIdeal.hpp>
#include <Reaktoro/Thermodynamics/Models/GaseousChemicalModelSpycherReed.hpp>
#include <Reaktoro/Thermodynamics/Models/GaseousChemicalModelSpycherPruessEnnis.hpp>

namespace Reaktoro {

struct GaseousPhase::Impl
{
    /// The gaseous mixture instance
    GaseousMixture mixture;

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
    setChemicalModelPengRobinson();
    setFluid();
}

GaseousPhase::~GaseousPhase()
{}

auto GaseousPhase::operator=(GaseousPhase other) -> GaseousPhase&
{
    Phase::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto GaseousPhase::setChemicalModelIdeal() -> GaseousPhase&
{
    PhaseChemicalModel model = gaseousChemicalModelIdeal(mixture());
    setChemicalModel(model);
    return *this;
}

auto GaseousPhase::setChemicalModelVanDerWaals() -> GaseousPhase&
{
    PhaseChemicalModel model = gaseousChemicalModelVanDerWaals(mixture());
    setChemicalModel(model);
    return *this;
}

auto GaseousPhase::setChemicalModelRedlichKwong() -> GaseousPhase&
{
    PhaseChemicalModel model = gaseousChemicalModelRedlichKwong(mixture());
    setChemicalModel(model);
    return *this;
}

auto GaseousPhase::setChemicalModelSoaveRedlichKwong() -> GaseousPhase&
{
    PhaseChemicalModel model = gaseousChemicalModelSoaveRedlichKwong(mixture());
    setChemicalModel(model);
    return *this;
}

auto GaseousPhase::setChemicalModelPengRobinson() -> GaseousPhase&
{
    PhaseChemicalModel model = gaseousChemicalModelPengRobinson(mixture());
    setChemicalModel(model);
    return *this;
}

auto GaseousPhase::setChemicalModelSpycherPruessEnnis() -> GaseousPhase&
{
    PhaseChemicalModel model = gaseousChemicalModelSpycherPruessEnnis(mixture());
    setChemicalModel(model);
    return *this;
}

auto GaseousPhase::setChemicalModelSpycherReed() -> GaseousPhase&
{
    PhaseChemicalModel model = gaseousChemicalModelSpycherReed(mixture());
    setChemicalModel(model);
    return *this;
}

auto GaseousPhase::mixture() const -> const GaseousMixture&
{
    return pimpl->mixture;
}

} // namespace Reaktoro
