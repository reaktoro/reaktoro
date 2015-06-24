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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Activity/MineralActivityIdeal.hpp>

namespace Reaktoro {
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
: MineralMixture(species), activity_fns(species.size())
{
    for(const auto& iter : species)
        setActivityModelIdeal(iter.name());
}

MineralPhase::MineralPhase(const MineralSpecies& species)
: MineralMixture(species), activity_fns(1)
{
    setActivityModelIdeal(species.name());
}

auto MineralPhase::setActivityModel(const std::string& species, const MineralActivityFunction& activity) -> void
{
    const Index ispecies = indexSpecies(species);
    if(ispecies < numSpecies())
        activity_fns[ispecies] = activity;
}

auto MineralPhase::setActivityModelIdeal(const std::string& species) -> void
{
    const Index ispecies = indexSpecies(species);

    if(ispecies < numSpecies())
        activity_fns[ispecies] = mineralActivityIdeal(species, *this);
}

auto MineralPhase::concentrations(double T, double P, const Vector& n) const -> ChemicalVector
{
    return molarFractions(n);
}

auto MineralPhase::activityConstants(double T, double P) const -> ThermoVector
{
    ThermoVector res(numSpecies());
    res.val.setConstant(1.0);
    return res;
}

auto MineralPhase::activityCoefficients(double T, double P, const Vector& n) const -> ChemicalVector
{
    return activities(T, P, n)/concentrations(T, P, n);
}

auto MineralPhase::activities(double T, double P, const Vector& n) const -> ChemicalVector
{
    MineralMixtureState mixture_state = state(T, P, n);
    const unsigned nspecies = numSpecies();
    ChemicalVector a(nspecies, nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        a.row(i) = activity_fns[i](mixture_state);
    return a;
}

} // namespace Reaktoro
