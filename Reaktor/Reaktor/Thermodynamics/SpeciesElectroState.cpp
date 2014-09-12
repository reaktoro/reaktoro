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

#include "SpeciesElectroState.hpp"

// Reaktor includes
#include <Reaktor/Species/AqueousSpecies.hpp>
#include <Reaktor/Thermodynamics/SpeciesElectroStateHKF.hpp>

namespace Reaktor {

SpeciesElectroState::SpeciesElectroState()
: reref(0), re(0), w(0), wT(0), wP(0), wTT(0), wTP(0), wPP(0)
{}

auto speciesElectro(double T, double P, const AqueousSpecies& species) -> SpeciesElectroState
{
    return speciesElectroHKF(T, P, species);
}

auto speciesElectro(const FunctionG& g, const AqueousSpecies& species) -> SpeciesElectroState
{
    return speciesElectroHKF(g, species);
}

} // namespace Reaktor
