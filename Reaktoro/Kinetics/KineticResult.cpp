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

#include <Reaktoro/Kinetics/KineticResult.hpp>

namespace Reaktoro {

auto KineticTiming::operator+=(const KineticTiming& other) -> KineticTiming&
{
    solve += other.solve;
    initialize += other.initialize;

    integrate += other.integrate;
    integrate_chemical_properties += other.integrate_chemical_properties;
    integrate_reaction_rates += other.integrate_reaction_rates;
    integrate_sensitivity += other.integrate_sensitivity;
    integrate_equilibration += other.integrate_equilibration;

    equilibrate += other.equilibrate;

    return *this;
}

auto KineticResult::operator+=(const KineticResult& other) -> KineticResult&
{
    timing += other.timing;
    equilibrium += other.equilibrium;
    smart_equilibrium += other.smart_equilibrium;
    return *this;
}
} // namespace Reaktoro