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
#include <ostream>
#include <map>
#include <tuple>
#include <vector>

// Reaktor includes
#include <Reaktor/Diagram/Utils/Types.hpp>
#include <Reaktor/Diagram/Utils/Segment.hpp>
#include <Reaktor/Diagram/Utils/Triangle.hpp>

namespace Reaktor {

class TernaryPlot
{
public:
	struct Options
	{
		/// The triangle length in refined areas of the mesh (default: 0.01)
		double fine_length;

		/// The triangle length in coarse areas of the mesh (default: 0.1)
		double coarse_length;

		/// The accuracy level for the calculation of a boundary point (default: 16)
		unsigned accuracy_level;

		/// The logical flag that indicates if the output of the ternary boundaries
		/// contains the indices of the points (default: true)
		bool output_indices;

		/// The logical flag that indicates if the output of the ternary boundaries
		/// contains the Barycentric coordinates of the points (default: true)
		bool output_barycentric;

		/// The logical flag that indicates if the output of the ternary boundaries
		/// contains the Cartesian coordinates of the points (default: true)
		bool output_cartesian;

		/// Constructs a default instance of @ref Options
		Options();
	};

	/**
	 * Constructs a default instance of @ref TernaryPlot
	 */
	TernaryPlot();

	/**
	 * Constructs an instance of @ref TernaryPlot
	 *
	 * @param func The ternary function for the ternary plot calculation
	 */
	TernaryPlot(const TernaryFunction& func);

	/**
	 * Sets the options for the calculations in the ternary plot
	 *
	 * @param options The options for the calculations
	 */
	auto setOptions(const Options& options) -> void;

	/**
	 * Gets the options for the ternary plot
	 */
	auto options() const -> const Options&;

	/**
	 * Gets all ternary boundaries found
	 */
	auto boundaries() const -> const TernaryBoundaries&;

	/**
	 * Builds the ternary plot by finding all ternary boundary lines
	 */
	auto build() -> void;

private:
	/// The ternary function
	TernaryFunction m_func;

	/// The ternary boundaries
	TernaryBoundaries m_boundaries;

	/// The options used for the calculations
	Options m_options;

	/// The ternary function converted to Cartesian coordinates
	std::function<std::string(const TernaryFunction&, const Point&)> m_funcpoint;

private:
	/**
	 * Recursively scans a boundary point in a segment
	 *
	 * The given segment is searched for a boundary point. The algorithm used is
	 * based on recursive bisection scheme for finding the root of nonlinear functions.
	 *
	 * If a boundary point is found, it is then inserted in @c boundaries.
	 *
	 * @param segment The segment where a boundary point is searched
	 * @param code_a The phase assemblage code at vertex @c a of @c segment
	 * @param code_b The phase assemblage code at vertex @c b of @c segment
	 */
	auto scanSegment(const Segment& segment, const std::string& code_a, const std::string& code_b) -> void;

	/**
	 * Recursively scans a boundary point in a triangle
	 *
	 * The given triangle is searched for several boundary points. The algorithm
	 * used is based on a recursive refinement scheme that yeilds smaller triangles
	 * along the boundary lines. At some point when the refined triangles satisfies
	 * a specifined fine mesh length, the algorithm @ref ScanSegment is used to seach
	 * for a boundary point along three internal segments, defined by the connection
	 * of the triangle vertices with its centroid.
	 *
	 * If boundary points are found, it is then inserted in @c boundaries.
	 *
	 * @param triangle The triangle where boundary points are searched
	 * @param code_a The phase assemblage code at vertex @c a of @c triangle
	 * @param code_b The phase assemblage code at vertex @c b of @c triangle
	 * @param code_c The phase assemblage code at vertex @c c of @c triangle
	 */
	auto scanTriangle(const Triangle& triangle, const std::string& code_a, const std::string& code_b, const std::string& code_c) -> void;

	/**
	 * Sorts the boundary points
	 *
	 * The boundary points in @c boundaries are sorted so that
	 * their order describes a smooth path, without overlappings.
	 */
	auto sort() -> void;
};

/**
 * Outputs the ternary plot
 *
 * @param out The output stream
 * @param plot The ternary plot as an instance of @ref TernaryPlot
 *
 * @return The updated output stream
 */
auto operator<<(std::ostream& out, const TernaryPlot& plot) -> std::ostream&;

}  /* namespace Reaktor */
