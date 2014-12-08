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
#include <memory>

// GEMS3K includes
#define IPMGEMPLUGIN
#define NOPARTICLEARRAY
#include <Reaktor/gems3k/node.h>

// Reaktor includes
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

/// An auxiliary type for a handle of the GEMS3K functionality
typedef std::shared_ptr<TNode> GemsHandle;

namespace internal {

/// Get the temperature of the GEMS3K handle (in units of K)
double temperature(const GemsHandle& node);

/// Get the pressure of the GEMS3K handle (in units of Pa)
double pressure(const GemsHandle& node);

unsigned numComponents(const GemsHandle& node);

unsigned numSpecies(const GemsHandle& node);

unsigned numPhases(const GemsHandle& node);

unsigned numSpeciesInPhase(const GemsHandle& node, unsigned iphase);

std::string componentName(const GemsHandle& node, unsigned icomponent);

std::string speciesName(const GemsHandle& node, unsigned ispecies);

std::string phaseName(const GemsHandle& node, unsigned iphase);

unsigned componentIndex(const GemsHandle& node, std::string component);

unsigned speciesIndex(const GemsHandle& node, std::string species);

unsigned phaseIndex(const GemsHandle& node, std::string phase);

double componentMolarMass(const GemsHandle& node, unsigned icomponent);

double speciesMolarMass(const GemsHandle& node, unsigned ispecies);

Vector speciesAmounts(const GemsHandle& node);

Vector componentAmounts(const GemsHandle& node);

Matrix formulaMatrix(const GemsHandle& node);

Vector standardChemicalPotentials(GemsHandle& node, double T, double P);

Vector chemicalPotentials(GemsHandle& node, double T, double P, const Vector& n);

bool converged(const GemsHandle& node);

} // namespace internal
} // namespace Reaktor
