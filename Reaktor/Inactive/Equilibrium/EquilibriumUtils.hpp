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

#pragma once

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {

// Reaktor forward declarations
class ChemicalState;
class Partitioning;

/**
 * Determines the existent equilibrium phases in the equilibrium state of a chemical system
 *
 * The existent phases at equilibrium are detected as follows. Let @c nphasei be the number of moles in the @a i-th phase.
 * Let @c ntotal be the total number of moles in the system. Let @c epsilon be a small scalar (i.e., 1.0e-6).
 * Then, if @c nphasei/@c ntotal > @c epsilon, we consider that the @a i-th phase is present at equilibrium.
 * @param state The state of the chemical system
 * @param partitioning The partitioning of the chemical system
 * @param epsilon The small scalar (default: 1.0e-6)
 * @return The indices of the existent equilibrium phases
 */
auto idxExistentPhases(const ChemicalState& state, const Partitioning& partitioning, double epsilon=1.0e-6) -> Indices;

} /* namespace Reaktor */
