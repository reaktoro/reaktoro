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

#include "Algorithms.hpp"

// C++ includes
#include <algorithm>
#include <cmath>

// Reaktor includes
#include <Reaktor/Diagram/Utils/Point.hpp>
#include <Reaktor/Diagram/Utils/Segment.hpp>

namespace Reaktor {

auto cross(const Segment& s1, const Segment& s2) -> bool
{
	const double& xa1 = s1.a.x;
	const double& xb1 = s1.b.x;
	const double& ya1 = s1.a.y;
	const double& yb1 = s1.b.y;
	const double& xa2 = s2.a.x;
	const double& xb2 = s2.b.x;
	const double& ya2 = s2.a.y;
	const double& yb2 = s2.b.y;

	const double num1 = (xb2 - xa2)*(ya2 - ya1) - (yb2 - ya2)*(xa2 - xa1);
	const double num2 = (xa1 - xa2)*(yb1 - ya1) - (ya1 - ya2)*(xb1 - xa1);
	const double den  = (xb2 - xa2)*(yb1 - ya1) - (xb1 - xa1)*(yb2 - ya2);

	if(std::abs(den) == 0.0) return false;

	const double t1 = num1/den;
	const double t2 = num2/den;

	return 0.0 < t1 and t1 < 1.0 and 0.0 < t2 and t2 < 1.0;
}

auto cross(const Segment& s, const std::vector<Point>& path) -> bool
{
	for(unsigned i = 1; i < path.size(); ++i)
		if(cross(s, Segment(path[i - 1], path[i])))
			return true;

	return false;
}

auto sort(const std::vector<Point>& points, const Point& ref) -> std::vector<Point>
{
    auto compare = [&](const Point& left, const Point& right)
    {
        const double dleft  = distanceSquared(ref, left);
        const double dright = distanceSquared(ref, right);

        return dleft < dright;
    };

    std::vector<Point> sorted(points);

    std::sort(sorted.begin(), sorted.end(), compare);

    return sorted;
}

auto overlapAtFront(const Point& point, const std::vector<Point>& path) -> bool
{
	return cross(Segment(path.front(), point), path);
}

auto overlapAtBack(const Point& point, const std::vector<Point>& path) -> bool
{
	return cross(Segment(path.back(), point), path);
}

auto pushFront(const Point& point, const std::vector<Point>& points) -> std::vector<Point>
{
	std::vector<Point> newpoints = { point };
    newpoints.insert(newpoints.end(), points.begin(), points.end());
    return newpoints;
}

auto pushBack(const Point& point, const std::vector<Point>& points) -> std::vector<Point>
{
	std::vector<Point> newpoints(points);
    newpoints.push_back(point);
    return newpoints;
}

auto deleteEntry(unsigned i, const std::vector<Point>& points) -> std::vector<Point>
{
	std::vector<Point> newpoints;
    newpoints.insert(newpoints.end(), points.begin(), points.begin() + i);
    newpoints.insert(newpoints.end(), points.begin() + i+1, points.end());
    return newpoints;
}

auto findPathFromFront(const std::vector<Point>& path, const std::vector<Point>& points, unsigned i) -> std::vector<Point>;

auto findPathFromBack(const std::vector<Point>& path, const std::vector<Point>& points, unsigned i) -> std::vector<Point>;

auto findPathHelper(const std::vector<Point>& path, const std::vector<Point>& points) -> std::vector<Point>
{
    if(points.empty())
        return path;

    std::vector<Point> sorted_from_front = sort(points, path.front());
    std::vector<Point> sorted_from_back  = sort(points, path.back());

    const double front_dist = distance(path.front(), sorted_from_front.front());
    const double back_dist  = distance(path.back(),  sorted_from_back.front());

    for(unsigned i = 0; i < points.size(); ++i)
    {
        if(front_dist < back_dist)
        {
            std::vector<Point> res_front = findPathFromFront(path, sorted_from_front, i);
            if(not res_front.empty()) return res_front;

            std::vector<Point> res_back = findPathFromBack(path, sorted_from_back, i);
            if(not res_back.empty()) return res_back;
        }
        else
        {
            std::vector<Point> res_back = findPathFromBack(path, sorted_from_back, i);
            if(not res_back.empty()) return res_back;

            std::vector<Point> res_front = findPathFromFront(path, sorted_from_front, i);
            if(not res_front.empty()) return res_front;
        }
    }

    return std::vector<Point>();
}

auto findPathFromFront(const std::vector<Point>& path, const std::vector<Point>& points, unsigned i) -> std::vector<Point>
{
    if(overlapAtFront(points[i], path))
        return std::vector<Point>();

    std::vector<Point> newpath = pushFront(points[i], path);
    std::vector<Point> newpoints = deleteEntry(i, points);

    return findPathHelper(newpath, newpoints);
}

auto findPathFromBack(const std::vector<Point>& path, const std::vector<Point>& points, unsigned i) -> std::vector<Point>
{
    if(overlapAtBack(points[i], path))
        return std::vector<Point>();

    std::vector<Point> newpath = pushBack(points[i], path);
    std::vector<Point> newpoints = deleteEntry(i, points);

    return findPathHelper(newpath, newpoints);
}

auto findPath(const std::vector<Point>& points) -> std::vector<Point>
{
    std::vector<Point> path = { points[0] };
    std::vector<Point> newpoints = deleteEntry(0, points);

    return findPathHelper(path, newpoints);
}

}  /* namespace Reaktor */
