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

#include "TernaryPlot.hpp"

// C++ includes
#include <iomanip>
#include <string>

// Reaktor includes
#include <Reaktor/Diagram/Utils/Algorithms.hpp>
#include <Reaktor/Diagram/Utils/BarycentricPoint.hpp>
#include <Reaktor/Diagram/Utils/Conversions.hpp>
#include <Reaktor/Diagram/Utils/Point.hpp>

namespace Reaktor {

TernaryPlot::Options::Options()
: fine_length(0.01), coarse_length(0.1),
  accuracy_level(16), output_indices(true),
  output_barycentric(true), output_cartesian(true)
{
}

TernaryPlot::TernaryPlot()
{
}

TernaryPlot::TernaryPlot(const TernaryFunction& func)
: m_func(func)
{
    m_funcpoint = [](const TernaryFunction& func, const Point& p)
    {
        BarycentricPoint q = barycentricPoint(p);
        return func(q.A, q.B, q.C);
    };
}

auto TernaryPlot::setOptions(const Options& options) -> void
{
    this->m_options = options;
}

auto TernaryPlot::options() const -> const Options&
{
    return m_options;
}

auto TernaryPlot::boundaries() const -> const TernaryBoundaries&
{
    return m_boundaries;
}

auto TernaryPlot::build() -> void
{
    m_boundaries.clear();

    const std::string code_a = m_func(1.0, 0.0, 0.0);
    const std::string code_b = m_func(0.0, 1.0, 0.0);
    const std::string code_c = m_func(0.0, 0.0, 1.0);

    const Point a = point(BarycentricPoint(1.0, 0.0, 0.0));
    const Point b = point(BarycentricPoint(0.0, 1.0, 0.0));
    const Point c = point(BarycentricPoint(0.0, 0.0, 1.0));

    scanTriangle(Triangle(a, b, c), code_a, code_b, code_c);

    sort();
}

auto TernaryPlot::scanSegment(const Segment& segment, const std::string& code_a, const std::string& code_b) -> void
{
    // Check if the phase assemblage codes at the vertices of the current segment are equal
    if(code_a == code_b)
        return;

    // Define an auxiliary lambda function that sorts the two phase assemblage codes
    auto sorted = [](const std::string& code1, const std::string& code2)
    {
        return std::make_tuple(std::min(code1, code2), std::max(code1, code2));
    };

    // The mid-point of the segment
    const Point mid = (segment.a + segment.b)/2.0;

    // Compute the length of the segments
    const double len = length(segment);

    // Check if the length of the current segment is less than the minimum scanning length
    if(len <= m_options.fine_length/m_options.accuracy_level)
    {
        m_boundaries[sorted(code_a, code_b)].push_back(barycentricPoint(mid));
    }
    else
    {
        // Evaluate the ternary function on the mid-point of the current segment
        const std::string code_mid = m_funcpoint(m_func, mid);

        scanSegment(Segment(segment.a, mid), code_a, code_mid);
        scanSegment(Segment(segment.b, mid), code_b, code_mid);
    }
}

auto TernaryPlot::scanTriangle(const Triangle& triangle, const std::string& code_a, const std::string& code_b, const std::string& code_c) -> void
{
    // Compute the length of the triangle
    const double len = length(triangle);

    // Check if the current triangle has a length inferior than the specified coarse mesh length
    // Also check if the phase assemblage codes at the vertices of the current triangle are equal
    if(len <= m_options.coarse_length) if(code_a == code_b and code_a == code_c and code_b == code_c)
            return;

    // Check if the current triangle has a length inferior that the specified fine mesh length
    if(len <= m_options.fine_length)
    {
        const Point d = (triangle.a + triangle.b + triangle.c)/3.0;

        const std::string code_d = m_funcpoint(m_func, d);

        const unsigned size = m_boundaries.size();

        scanSegment(Segment(triangle.a, d), code_a, code_d); if(size != m_boundaries.size()) return;
        scanSegment(Segment(triangle.b, d), code_b, code_d); if(size != m_boundaries.size()) return;
        scanSegment(Segment(triangle.c, d), code_c, code_d); if(size != m_boundaries.size()) return;
    }
    else
    {
        // Calculate the mid-points along the edges `ab`, `bc` and `ca`
        Point mab = (triangle.a + triangle.b)/2.0;
        Point mbc = (triangle.b + triangle.c)/2.0;
        Point mca = (triangle.c + triangle.a)/2.0;

        // Calculate the phase assemblage codes at the mid-points of the triangle
        const std::string code_mab = m_funcpoint(m_func, mab);
        const std::string code_mbc = m_funcpoint(m_func, mbc);
        const std::string code_mca = m_funcpoint(m_func, mca);

        // Refine the scan for boundary points in the four refined triangles
        scanTriangle(Triangle(triangle.a, mab, mca), code_a, code_mab, code_mca);
        scanTriangle(Triangle(triangle.b, mbc, mab), code_b, code_mbc, code_mab);
        scanTriangle(Triangle(triangle.c, mca, mbc), code_c, code_mca, code_mbc);
        scanTriangle(Triangle(mab, mbc, mca), code_mab, code_mbc, code_mca);
    }
}

auto sortHelper(const TernaryBoundary& boundary) -> TernaryBoundary
{
    std::vector<Point> points;
    points.reserve(boundary.size());
    for(const BarycentricPoint& q : boundary)
        points.push_back(point(q));

    points = findPath(points);

    TernaryBoundary new_boundary;
    new_boundary.reserve(boundary.size());
    for(const Point& p : points)
        new_boundary.push_back(barycentricPoint(p));

    return new_boundary;
}

auto TernaryPlot::sort() -> void
{
    for(const auto& pair : m_boundaries)
        m_boundaries[pair.first] = sortHelper(pair.second);
}

auto operator<<(std::ostream& out, const TernaryPlot& plot) -> std::ostream&
{
    const TernaryPlot::Options& options = plot.options();
    const TernaryBoundaries& boundaries = plot.boundaries();

    unsigned nfill = 0;
    nfill += options.output_barycentric ? 45 : 0;
    nfill += options.output_cartesian   ? 30 : 0;
    nfill += options.output_indices     ?  5 : 0;

    const std::string bar(nfill, '=');

    out << bar << std::endl;
    if(options.output_barycentric) out << std::setw(15) << std::left << "A";
    if(options.output_barycentric) out << std::setw(15) << std::left << "B";
    if(options.output_barycentric) out << std::setw(15) << std::left << "C";
    if(options.output_cartesian)   out << std::setw(15) << std::left << "x";
    if(options.output_cartesian)   out << std::setw(15) << std::left << "y";
    if(options.output_indices)     out << std::setw(5)  << std::left << "i";
    out << std::endl << bar << std::endl;

    for(const auto& pair : boundaries)
    {
        std::string code1, code2; std::tie(code1, code2) = pair.first;
        const TernaryBoundary& boundary = pair.second;

        out << "PhaseTransition: " << code1 << "-" << code2 << std::endl;
        out << bar << std::endl;

        for(unsigned i = 0; i < boundary.size(); ++i)
        {
            const Point pnt = point(boundary[i]);

            if(options.output_barycentric) out << std::setw(15) << std::left << boundary[i].A;
            if(options.output_barycentric) out << std::setw(15) << std::left << boundary[i].B;
            if(options.output_barycentric) out << std::setw(15) << std::left << boundary[i].C;
            if(options.output_cartesian)   out << std::setw(15) << std::left << pnt.x;
            if(options.output_cartesian)   out << std::setw(15) << std::left << pnt.y;
            if(options.output_indices)     out << std::setw(5)  << std::left << i;
            out << std::endl;
        }

        out << std::endl;
    }

    return out;
}

}  /* namespace Reaktor */
