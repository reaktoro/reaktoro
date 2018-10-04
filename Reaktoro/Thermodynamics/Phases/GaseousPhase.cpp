// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

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
: Phase(), pimpl(new Impl())
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
    setType(PhaseType::Gas);
    setSpecies(species);
    setChemicalModelPengRobinson();
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
