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

#pragma once

// C++ includes
#include <functional>

// Reaktor includes
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class ScalarResult;
class VectorResult;

/// Define the function signature of the standart chemical potential of a species (in units of J/mol)
/// @param T The temperature value (in units of K)
/// @param P The pressure value (in units of Pa)
/// @return The standart chemical potential of the species (in units of J/mol)
/// @see Species
/// @ingroup Core
typedef std::function<
    double(double T, double P)>
        ChemicalPotential;

/// Define the function signature of the equilibrium constant of a reaction
/// @param T The temperature value (in units of K)
/// @param P The pressure value (in units of Pa)
/// @return The equilibrium constant of the reaction
/// @see Reaction
/// @ingroup Core
typedef std::function<
    double(double T, double P)>
        EquilibriumConstant;

/// Define the function signature of the concentration of a phase
/// @param n The molar amounts of the species in the phase (in units of mol)
/// @return The concentrations of the species in the phase
/// @see Phase
/// @ingroup Core
typedef std::function<
    Vector(const SubVector& n)>
        Concentration;

/// Define the function signature of the activity of a species in a phase
/// @param T The temperature value (in units of K)
/// @param P The pressure value (in units of Pa)
/// @param n The molar amounts of the species in the phase (in units of mol)
/// @return The activities of the species in the phase and their molar derivatives
/// @see Phase
/// @ingroup Core
typedef std::function<
    ScalarResult(double T, double P, const SubVector& n)>
        Activity;

/// Define the function signature of the density of the phase (in units of kg/m3)
/// @param T The temperature value (in units of K)
/// @param P The pressure value (in units of Pa)
/// @param n The molar amounts of the species in the phase (in units of mol)
/// @return The density of the phase and its molar derivatives (in units of kg/m3)
/// @ingroup Core
typedef std::function<
    ScalarResult(double T, double P, const SubVector& n)>
        Density;

/// Define the function signature of the rate of a reaction (in units of mol/s)
/// @param T The temperature value (in units of K)
/// @param P The pressure value (in units of Pa)
/// @param n The molar amounts of all species (in units of mol)
/// @param a The activities of all species and their molar derivatives
/// @return The rate of the reaction and its molar derivatives (in units of mol/s)
/// @see Reaction
/// @ingroup Core
typedef std::function<
    ScalarResult(double T, double P, const Vector& n, const VectorResult& a)>
        Rate;

} /* namespace Reaktor */
