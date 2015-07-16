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

#include "MineralPhase.hpp"

// C++ includes
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Activity/MineralActivityIdeal.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/MineralMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/MineralChemicalModelIdeal.hpp>

namespace Reaktoro {
namespace internal {

auto nameMineralPhase(const MineralMixture& mixture) -> std::string
{
    std::stringstream name;
    for(const MineralSpecies& iter : mixture.species())
        name << iter.name() << "-";
    std::string str = name.str();
    str = str.substr(0, str.size() - 1);
    return str;
}

} // namespace internal

struct MineralPhase::Impl
{
    /// The mineral mixture instance
    MineralMixture mixture;

    /// The functions that calculates the activities of selected species
    std::map<Index, MineralActivityFunction> activities;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance
    Impl(const MineralMixture& mixture)
    : mixture(mixture)
    {}
};

MineralPhase::MineralPhase()
: Phase()
{}

MineralPhase::MineralPhase(const MineralPhase& other)
: Phase(other), pimpl(new Impl(*other.pimpl))
{}

MineralPhase::MineralPhase(const MineralMixture& mixture)
: pimpl(new Impl(mixture))
{
    // Convert the MineralSpecies instances to Species instances
    std::vector<Species> species;
    for(const MineralSpecies& x : mixture.species())
        species.push_back(x);

    // Set the Phase attributes
    setName(internal::nameMineralPhase(mixture));
    setSpecies(species);
    setReferenceState(PhaseReferenceState::IdealSolution);
    setChemicalModelIdeal();
}

MineralPhase::MineralPhase(const MineralSpecies& species)
: MineralPhase(MineralMixture(std::vector<MineralSpecies>{species}))
{}

MineralPhase::~MineralPhase()
{}

auto MineralPhase::operator=(MineralPhase other) -> MineralPhase&
{
    Phase::operator=(other);
    pimpl = std::move(other.pimpl);
    return *this;
}

auto MineralPhase::setChemicalModelIdeal() -> void
{
    // Create the aqueous chemical model
    PhaseChemicalModel model = mineralChemicalModelIdeal(mixture());

    setChemicalModel(model);
}

auto MineralPhase::setActivityModel(const std::string& species, const MineralActivityFunction& activity) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = activity;
}

auto MineralPhase::setActivityModelIdeal(const std::string& species) -> void
{
    const Index ispecies = indexSpecies(species);

    if(ispecies < numSpecies())
        pimpl->activities[ispecies] = mineralActivityIdeal(species, mixture());
}

auto MineralPhase::mixture() const -> const MineralMixture&
{
    return pimpl->mixture;
}

} // namespace Reaktoro
