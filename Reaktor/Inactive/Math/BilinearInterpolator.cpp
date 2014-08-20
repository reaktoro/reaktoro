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

#include "BilinearInterpolator.hpp"

// C++ includes
#include <cmath>
#include <iomanip>
#include <sstream>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Math/BilinearInterpolator.hpp>

namespace Reaktor {
namespace internal {

auto binarySearchHelper(double p, const std::vector<double>& coordinates, unsigned begin, unsigned end) -> unsigned
{
    if(end - begin == 1)
        return begin;

    unsigned mid = (begin + end)/2;

    if(p < coordinates[mid])
        return binarySearchHelper(p, coordinates, begin, mid);
    else
        return binarySearchHelper(p, coordinates, mid, end);
}

auto binarySearch(double p, const std::vector<double>& coordinates) -> unsigned
{
    return binarySearchHelper(p, coordinates, 0, coordinates.size());
}

inline auto interpolationOutOfBoundsError(double x, double xA, double xB, double y, double yA, double yB) -> void
{
    Exception exception;
    exception.error << "Unable to perform an interpolation at the coordinate pair (" << x << ", " << y << ").";
    exception.reason << "Either the x- or y-coordinate is out of bound, where " << xA << " < x < " << xB << " and " << yA << " < y < " << yB << ".";
    raise(exception);
}

} /* namespace internal */

BilinearInterpolator::BilinearInterpolator()
{}

BilinearInterpolator::BilinearInterpolator(
	const std::vector<double>& xcoordinates,
	const std::vector<double>& ycoordinates,
	const std::vector<double>& data)
: xcoordinates$(xcoordinates),
  ycoordinates$(ycoordinates),
  data$(data)
{}

BilinearInterpolator::BilinearInterpolator(
    const std::vector<double>& xcoordinates,
    const std::vector<double>& ycoordinates,
    const std::function<double(double, double)>& function)
: xcoordinates$(xcoordinates),
  ycoordinates$(ycoordinates),
  data$(xcoordinates.size() * ycoordinates.size())
{
    unsigned k = 0;
    for(unsigned j = 0; j < ycoordinates.size(); ++j)
        for(unsigned i = 0; i < xcoordinates.size(); ++i, ++k)
            data$[k] = function(xcoordinates[i], ycoordinates[j]);
}

auto BilinearInterpolator::setCoordinatesX(const std::vector<double>& xcoordinates) -> void
{
    xcoordinates$ = xcoordinates;
}

auto BilinearInterpolator::setCoordinatesY(const std::vector<double>& ycoordinates) -> void
{
    ycoordinates$ = ycoordinates;
}

auto BilinearInterpolator::setData(const std::vector<double>& data) -> void
{
    data$ = data;
}

auto BilinearInterpolator::xCoodinates() const -> const std::vector<double>&
{
	return xcoordinates$;
}

auto BilinearInterpolator::yCoodinates() const -> const std::vector<double>&
{
	return ycoordinates$;
}

auto BilinearInterpolator::data() const -> const std::vector<double>&
{
	return data$;
}

auto BilinearInterpolator::initialised() const -> bool
{
    return not data$.empty();
}

auto BilinearInterpolator::operator()(double x, double y) const -> double
{
    const double xA = xcoordinates$.front();
    const double xB = xcoordinates$.back();
    const double yA = ycoordinates$.front();
    const double yB = ycoordinates$.back();

    x = std::max(xA, std::min(x, xB));
    y = std::max(yA, std::min(y, yB));

	const unsigned sizex = xcoordinates$.size();
	const unsigned sizey = ycoordinates$.size();

	const double i = internal::binarySearch(x, xcoordinates$);
	const double j = internal::binarySearch(y, ycoordinates$);

	const auto k = [=](unsigned i, unsigned j) { return i + j*sizex; };

	if(i == sizex or j == sizey)
	    internal::interpolationOutOfBoundsError(x, xA, xB, y, yA, yB);

	const double x1 = xcoordinates$[i];
	const double x2 = xcoordinates$[i + 1];

	const double y1 = ycoordinates$[j];
	const double y2 = ycoordinates$[j + 1];

	const double z11 = data$[k(i  , j  )]; // z at (x1, y1)
    const double z21 = data$[k(i+1, j  )]; // z at (x2, y1)
    const double z12 = data$[k(i  , j+1)]; // z at (x1, y2)
    const double z22 = data$[k(i+1, j+1)]; // z at (x2, y2)

	const double f11 =  z11*(x2 - x)*(y2 - y);
	const double f12 = -z12*(x2 - x)*(y1 - y);
	const double f21 = -z21*(x1 - x)*(y2 - y);
	const double f22 =  z22*(x1 - x)*(y1 - y);

	return (f11 + f12 + f21 + f22)/((x2 - x1)*(y2 - y1));
}

auto operator<<(std::ostream& out, const BilinearInterpolator& interpolator) -> std::ostream&
{
	const auto& xcoordinates = interpolator.xCoodinates();
	const auto& ycoordinates = interpolator.yCoodinates();
	const auto& data         = interpolator.data();

	const unsigned sizex = xcoordinates.size();
    const unsigned sizey = ycoordinates.size();

	out << std::setw(15) << std::right << "y/x";
	for(auto x : xcoordinates)
		out << std::setw(15) << std::right << x;
	out << std::endl;

	for(unsigned j = 0; j < sizey; ++j)
	{
		out << std::setw(15) << std::right << ycoordinates[j];
		for(unsigned i = 0; i < sizex; ++i)
			out << std::setw(15) << std::right << data[i + j*sizex];
		out << std::endl;
	}

	return out;
}

} /* namespace Reaktor */
