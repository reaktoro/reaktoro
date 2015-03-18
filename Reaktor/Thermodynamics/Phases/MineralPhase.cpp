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

#include "MineralPhase.hpp"

// C++ includes
#include <sstream>

// Reaktor includes
#include <Reaktor/Activity/MineralActivityIdeal.hpp>
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Core/Phase.hpp>

namespace Reaktor {
namespace internal {

auto nameMineralPhase(const MineralPhase& phase) -> std::string
{
    std::stringstream name;
    for(const MineralSpecies& iter : phase.species())
        name << iter.name() << "-";
    std::string str = name.str();
    str = str.substr(0, str.size() - 1);
    return str;
}

} /* namespace internal */

MineralPhase::MineralPhase()
: MineralMixture()
{}

MineralPhase::MineralPhase(const std::vector<MineralSpecies>& species)
: MineralMixture(species), activities$(species.size())
{
    for(const auto& iter : species)
        setActivityModelIdeal(iter.name());
}

auto MineralPhase::setActivityModel(const std::string& species, const MineralActivity& activity) -> void
{
    const Index ispecies = idxSpecies(species);
    if(ispecies < numSpecies())
        activities$[ispecies] = activity;
}

auto MineralPhase::setActivityModelIdeal(const std::string& species) -> void
{
    const Index ispecies = idxSpecies(species);

    if(ispecies < numSpecies())
        activities$[ispecies] = mineralActivityIdeal(species, *this);
}

auto MineralPhase::params(double T, double P, const Vector& n) const -> MineralActivityParams
{
    MineralActivityParams params;

    params.T = T;
    params.P = P;
    params.n = n;
    params.x = molarFractions(n);

    return params;
}

auto MineralPhase::concentrations(const Vector& n) const -> Vector
{
    // The total amount of moles in the mineral phase
    const double ntotal = n.sum();

    // Check if the phase has zero number of moles
    if(ntotal == 0.0) return zeros(n.rows());

    // Calculate the molar fractions of the mineral species
    return n/ntotal;
}

auto MineralPhase::activities(double T, double P, const Vector& n) const -> PartialVector
{
    MineralActivityParams pars = params(T, P, n);

    const unsigned N = numSpecies();

    PartialVector a = partialVector(zeros(N), zeros(N, N));

    for(unsigned i = 0; i < N; ++i)
    {
        const PartialScalar res = activities$[i](pars);

        func(a)[i] = func(res);
        grad(a).row(i) = grad(res);
    }

    return a;
}

auto createPhase(const MineralPhase& phase) -> Phase
{
    // Create the mineral species as Species instances
    std::vector<Species> species;
    for(const MineralSpecies& iter : phase.species())
        species.push_back(createSpecies(iter));

    // Define the concentration function of the mineral phase
    Concentration concentration = [=](const Vector& n) -> Vector
    {
        return phase.concentrations(n);
    };

    // Define the activity function of the mineral phase
    Activity activity = [=](double T, double P, const Vector& n)
    {
        return phase.activities(T, P, n);
    };

    Phase converted;
    converted.setName(internal::nameMineralPhase(phase));
    converted.setSpecies(species);
    converted.setConcentration(concentration);
    converted.setActivity(activity);

    return converted;
}

} // namespace Reaktor
