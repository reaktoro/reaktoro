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

// C++ includes
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace Reaktor {

// Reaktor forward declarations
struct BarycentricPoint;
struct Point;

/**
 * Typedef that defines the phase transition across a phase boundary
 *
 * The phase transition across a phase boundary is defined as a pair of
 * codes (code1, code2) that indicates the change in the phase assemblage
 * state across a phase boundary.
 */
using PhaseTransition = std::tuple<std::string, std::string>;

/**
 * Typedef that defines the signature of phase functions
 *
 * A phase function f = f(x, y) is a function that returns the
 * phase assemblage code at the Cartesian point (x, y).
 */
using PhaseFunction = std::function<std::string(double, double)>;

/**
 * Typedef that defines a phase boundary
 *
 * A phase boundary is a container of Cartesian points (x, y)
 * describing a path across which a phase state transition is observed.
 */
using PhaseBoundary = std::vector<Point>;

/**
 * Typedef that defines a set of phase boundaries
 *
 * The phase boundaries are arranged as a map so that they are accessed
 * through the use of a key representing the phase state transition across
 * the boundary.
 */
using PhaseBoundaries = std::map<PhaseTransition, PhaseBoundary>;

/**
 * Typedef that defines the signature of ternary functions
 *
 * A ternary phase function f = f(A, B, C) is a function that
 * returns the phase assemblage code at the barycentric point
 * (A, B, C), where A + B + C = 1 and 0 <= A, B, C <= 1.
 */
using TernaryFunction = std::function<std::string(double, double, double)>;

/**
 * Typedef that defines a ternary phase boundary
 *
 * A ternary boundary is a container of barycentric points (A, B, C)
 * describing a path across which a phase state transition is observed.
 */
using TernaryBoundary = std::vector<BarycentricPoint>;

/**
 * Typedef that defines a set of ternary boundaries
 *
 * The ternary boundaries are arranged as a map so that they are accessed
 * through the use of a key representing the phase state transition across
 * the boundary.
 */
using TernaryBoundaries = std::map<PhaseTransition, TernaryBoundary>;

}  /* namespace Reaktor */
