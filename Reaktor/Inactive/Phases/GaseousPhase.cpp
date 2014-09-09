/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "GaseousPhase.hpp"

// Reaktor includes
#include <Reaktor/Activity/GaseousActivityDuanSun.hpp>
#include <Reaktor/Activity/GaseousActivityIdeal.hpp>
#include <Reaktor/Activity/GaseousActivityPengRobinson.hpp>
#include <Reaktor/Activity/GaseousActivitySpycherPruess.hpp>
#include <Reaktor/Activity/GaseousActivitySpycherReed.hpp>
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

GaseousPhase::GaseousPhase()
: GaseousMixture()
{}

GaseousPhase::GaseousPhase(const std::vector<GaseousSpecies>& species)
: GaseousMixture(species), activities$(species.size())
{
    for(const auto& iter : species)
        setActivityModelIdeal(iter.name());
}

auto GaseousPhase::setActivityModel(const std::string& species, const GaseousActivity& activity) -> void
{
    const Index ispecies = idxSpecies(species);
    if(ispecies < numSpecies())
        activities$[ispecies] = activity;
}

auto GaseousPhase::setActivityModelIdeal(const std::string& species) -> void
{
    const Index ispecies = idxSpecies(species);

    if(ispecies < numSpecies())
        activities$[ispecies] = gaseousActivityIdeal(species, *this);
}

auto GaseousPhase::setActivityModelDuanSunCO2() -> void
{
    const Index ispecies = idxSpecies("CO2(g)");

    if(ispecies < numSpecies())
        activities$[ispecies] = gaseousActivityDuanSunCO2(*this);
}

auto GaseousPhase::setActivityModelSpycherPruessH2OCO2() -> void
{
    const Index iH2O = idxSpecies("H2O(g)");
    const Index iCO2 = idxSpecies("CO2(g)");

    const auto& functions = gaseousActivitySpycherPruessH2OCO2(*this);

    if(iH2O < numSpecies()) activities$[iH2O] = functions[0];
    if(iCO2 < numSpecies()) activities$[iCO2] = functions[1];
}

auto GaseousPhase::setActivityModelSpycherReedH2OCO2CH4() -> void
{
    const Index iH2O = idxSpecies("H2O(g)");
    const Index iCO2 = idxSpecies("CO2(g)");
    const Index iCH4 = idxSpecies("CH4(g)");

    const auto& functions = gaseousActivitySpycherReedH2OCO2CH4(*this);

    if(iH2O < numSpecies()) activities$[iH2O] = functions[0];
    if(iCO2 < numSpecies()) activities$[iCO2] = functions[1];
    if(iCH4 < numSpecies()) activities$[iCH4] = functions[2];
}

auto GaseousPhase::setActivityModelPengRobinson(const std::string& species) -> void
{
    const Index idx_species = idxSpecies(species);

    if(idx_species < numSpecies())
        activities$[idx_species] = gaseousActivityPengRobinson(species, *this);
}

auto GaseousPhase::params(double T, double P, const Vector& n) const -> GaseousActivityParams
{
    GaseousActivityParams params;

    params.T = T;
    params.P = P;
    params.n = n;
    params.x = molarFractions(n);

    return params;
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

auto GaseousPhase::activities(double T, double P, const Vector& n) const -> VectorResult
{
    GaseousActivityParams pars = params(T, P, n);

    const unsigned N = numSpecies();

    VectorResult a = partialVector(zeros(N), zeros(N, N));

    for(unsigned i = 0; i < N; ++i)
    {
        const ScalarResult res = activities$[i](pars);

        func(a)[i] = func(res);
        grad(a).row(i) = grad(res);
    }

    return a;
}

} // namespace Reaktor
