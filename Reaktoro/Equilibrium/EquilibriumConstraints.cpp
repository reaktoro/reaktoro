// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "EquilibriumConstraints.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

EquilibriumConstraints::EquilibriumConstraints()
{}

auto EquilibriumConstraints::temperature(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::temperature has not been implemented yet.");
}

auto EquilibriumConstraints::pressure(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::pressure has not been implemented yet.");
}

auto EquilibriumConstraints::volume(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::volume has not been implemented yet.");
}

auto EquilibriumConstraints::internalEnergy(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::internalEnergy has not been implemented yet.");
}

auto EquilibriumConstraints::enthalpy(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::enthalpy has not been implemented yet.");
}

auto EquilibriumConstraints::gibbsEnergy(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::gibbsEnergy has not been implemented yet.");
}

auto EquilibriumConstraints::helmholtzEnergy(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::helmholtzEnergy has not been implemented yet.");
}

auto EquilibriumConstraints::entropy(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::entropy has not been implemented yet.");
}

auto EquilibriumConstraints::chemicalPotential(String species, real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::chemicalPotential has not been implemented yet.");
}

auto EquilibriumConstraints::chemicalPotential(String species, Fn<real(real,real)> fn, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::chemicalPotential has not been implemented yet.");
}

auto EquilibriumConstraints::activity(String species, real value) -> void
{
    error(true, "Method EquilibriumConstraints::activity has not been implemented yet.");
}

auto EquilibriumConstraints::fugacity(String species, real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::fugacity has not been implemented yet.");
}

auto EquilibriumConstraints::pH(real value) -> void
{
    error(true, "Method EquilibriumConstraints::pH has not been implemented yet.");
}

auto EquilibriumConstraints::pe(real value) -> void
{
    error(true, "Method EquilibriumConstraints::pe has not been implemented yet.");
}

auto EquilibriumConstraints::Eh(real value, String unit) -> void
{
    error(true, "Method EquilibriumConstraints::Eh has not been implemented yet.");
}

} // namespace Reaktoro
