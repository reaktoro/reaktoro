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

#include "BilinearInterpolator.hpp"

// C++ includes
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>

namespace Reaktoro {
namespace {

auto binarySearchHelper(double p, const std::vector<double>& coordinates, size_t begin, size_t end) -> size_t
{
    if(end - begin == 1)
        return begin;

    auto mid = (begin + end)/2;

    if(p < coordinates[mid])
        return binarySearchHelper(p, coordinates, begin, mid);
    else
        return binarySearchHelper(p, coordinates, mid, end);
}

auto binarySearch(double p, const std::vector<double>& coordinates) -> size_t
{
    return binarySearchHelper(p, coordinates, 0, coordinates.size());
}

} // namespace

BilinearInterpolator::BilinearInterpolator()
{}

BilinearInterpolator::BilinearInterpolator(
    const std::vector<double>& xcoordinates,
    const std::vector<double>& ycoordinates,
    const std::vector<double>& data)
: m_xcoordinates(xcoordinates),
  m_ycoordinates(ycoordinates),
  m_data(data)
{}

BilinearInterpolator::BilinearInterpolator(
    const std::vector<double>& xcoordinates,
    const std::vector<double>& ycoordinates,
    const std::function<double(double, double)>& function)
: m_xcoordinates(xcoordinates),
  m_ycoordinates(ycoordinates),
  m_data(xcoordinates.size() * ycoordinates.size())
{
    unsigned k = 0;
    for(unsigned j = 0; j < ycoordinates.size(); ++j)
        for(unsigned i = 0; i < xcoordinates.size(); ++i, ++k)
            m_data[k] = function(xcoordinates[i], ycoordinates[j]);
}

auto BilinearInterpolator::setCoordinatesX(const std::vector<double>& xcoordinates) -> void
{
    m_xcoordinates = xcoordinates;
}

auto BilinearInterpolator::setCoordinatesY(const std::vector<double>& ycoordinates) -> void
{
    m_ycoordinates = ycoordinates;
}

auto BilinearInterpolator::setData(const std::vector<double>& data) -> void
{
    m_data = data;
}

auto BilinearInterpolator::xCoordinates() const -> const std::vector<double>&
{
    return m_xcoordinates;
}

auto BilinearInterpolator::yCoordinates() const -> const std::vector<double>&
{
    return m_ycoordinates;
}

auto BilinearInterpolator::data() const -> const std::vector<double>&
{
    return m_data;
}

auto BilinearInterpolator::empty() const -> bool
{
    return m_data.empty();
}

auto BilinearInterpolator::operator()(double x, double y) const -> double
{
    // Check if the interpolation data contains only one point
    if(m_data.size() == 1) return m_data[0];

    const double xA = m_xcoordinates.front();
    const double xB = m_xcoordinates.back();
    const double yA = m_ycoordinates.front();
    const double yB = m_ycoordinates.back();

    x = std::max(xA, std::min(x, xB));
    y = std::max(yA, std::min(y, yB));

    const auto size_x = m_xcoordinates.size();
    const auto size_y = m_ycoordinates.size();

    const auto index_x = binarySearch(x, m_xcoordinates);
    const auto i = index_x == size_x - 1 ? index_x - 1 : index_x;

    const auto index_y = binarySearch(y, m_ycoordinates);
    const auto j = index_y == size_y - 1 ? index_y - 1 : index_y;

    const auto k = [=](size_t i, size_t j) { return i + j*size_x; };

    const double x1 = m_xcoordinates[i];
    const double x2 = m_xcoordinates[i + 1];

    const double y1 = m_ycoordinates[j];
    const double y2 = m_ycoordinates[j + 1];

    const double z11 = m_data[k(i  , j  )]; // z at (x1, y1)
    const double z21 = m_data[k(i+1, j  )]; // z at (x2, y1)
    const double z12 = m_data[k(i  , j+1)]; // z at (x1, y2)
    const double z22 = m_data[k(i+1, j+1)]; // z at (x2, y2)

    const double f11 =  z11*(x2 - x)*(y2 - y);
    const double f12 = -z12*(x2 - x)*(y1 - y);
    const double f21 = -z21*(x1 - x)*(y2 - y);
    const double f22 =  z22*(x1 - x)*(y1 - y);

    return (f11 + f12 + f21 + f22)/((x2 - x1)*(y2 - y1));
}

auto operator<<(std::ostream& out, const BilinearInterpolator& interpolator) -> std::ostream&
{
    const auto& xcoordinates = interpolator.xCoordinates();
    const auto& ycoordinates = interpolator.yCoordinates();
    const auto& data         = interpolator.data();

    const auto sizex = xcoordinates.size();
    const auto sizey = ycoordinates.size();

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

} // namespace Reaktoro
