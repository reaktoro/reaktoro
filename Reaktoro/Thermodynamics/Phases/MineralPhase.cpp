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

#include "MineralPhase.hpp"

// C++ includes
#include <string>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/MineralMixture.hpp>
#include <Reaktoro/Thermodynamics/Models/MineralChemicalModelIdeal.hpp>
#include <Reaktoro/Thermodynamics/Models/MineralChemicalModelRedlichKister.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

namespace Reaktoro {
namespace internal {

auto nameMineralPhase(const MineralMixture& mixture) -> std::string
{
    std::string name;
    for(const MineralSpecies& iter : mixture.species())
        name += iter.name() + "-";
    name = name.substr(0, name.size() - 1);
    return name;
}

} // namespace internal

struct MineralPhase::Impl
{
    /// The mineral mixture instance
    MineralMixture mixture;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a custom Impl instance
    Impl(const MineralMixture& mixture)
    : mixture(mixture)
    {}
};

MineralPhase::MineralPhase()
: Phase(), pimpl(new Impl())
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
    setType(PhaseType::Solid);
    setSpecies(species);
    setChemicalModelIdeal();
}

MineralPhase::MineralPhase(const MineralSpecies& species)
: MineralPhase(MineralMixture(std::vector<MineralSpecies>{species}))
{}

auto MineralPhase::setChemicalModelIdeal() -> MineralPhase&
{
    PhaseChemicalModel model = mineralChemicalModelIdeal(mixture());
    setChemicalModel(model);
    return *this;
}

auto MineralPhase::setChemicalModelRedlichKister(double a0, double a1, double a2) -> MineralPhase&
{
    PhaseChemicalModel model = mineralChemicalModelRedlichKister(mixture(), a0, a1, a2);
    setChemicalModel(model);
    return *this;
}

auto MineralPhase::mixture() const -> const MineralMixture&
{
    return pimpl->mixture;
}

} // namespace Reaktoro
