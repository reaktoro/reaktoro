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
    PhaseChemicalModel model = gaseousChemicalModelIdeal(mixture());
    setChemicalModel(model);
}

auto GaseousPhase::setChemicalModelVanDerWaals() -> void
{
    PhaseChemicalModel model = gaseousChemicalModelVanDerWaals(mixture());
    setChemicalModel(model);
}

auto GaseousPhase::setChemicalModelRedlichKwong() -> void
{
    PhaseChemicalModel model = gaseousChemicalModelRedlichKwong(mixture());
    setChemicalModel(model);
}

auto GaseousPhase::setChemicalModelSoaveRedlichKwong() -> void
{
    PhaseChemicalModel model = gaseousChemicalModelSoaveRedlichKwong(mixture());
    setChemicalModel(model);
}

auto GaseousPhase::setChemicalModelPengRobinson() -> void
{
    PhaseChemicalModel model = gaseousChemicalModelPengRobinson(mixture());
    setChemicalModel(model);
}

auto GaseousPhase::setChemicalModelSpycherPruessEnnis() -> void
{
    PhaseChemicalModel model = gaseousChemicalModelSpycherPruessEnnis(mixture());
    setChemicalModel(model);
}

auto GaseousPhase::setChemicalModelSpycherReed() -> void
{
    PhaseChemicalModel model = gaseousChemicalModelSpycherReed(mixture());
    setChemicalModel(model);
}

auto GaseousPhase::mixture() const -> const GaseousMixture&
{
    return pimpl->mixture;
}

} // namespace Reaktoro
