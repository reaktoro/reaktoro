// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

/// A class used to calculate bilinear interpolation of data in two dimensions.
class BilinearInterpolator
{
public:
    /// Construct a default BilinearInterpolator instance
    BilinearInterpolator();

    /// Construct a BilinearInterpolator instance with given data
    /// @param xcoordinates The x-coordinates for the interpolation
    /// @param ycoordinates The y-coordinates for the interpolation
    /// @param data The data to be interpolated over the (x, y) coordinates
    BilinearInterpolator(
        const Vec<double>& xcoordinates,
        const Vec<double>& ycoordinates,
        const Vec<double>& data);

    /// Construct a BilinearInterpolator instance with given data
    /// @param xcoordinates The x-coordinates for the interpolation
    /// @param ycoordinates The y-coordinates for the interpolation
    /// @param data The data to be interpolated over the (x, y) coordinates
    BilinearInterpolator(
        const Vec<double>& xcoordinates,
        const Vec<double>& ycoordinates,
        const Vec<Vec<double>>& data);

    /// Construct a BilinearInterpolator instance with given function
    /// @param xcoordinates The x-coordinates for the interpolation
    /// @param ycoordinates The y-coordinates for the interpolation
    /// @param function The function to be interpolated over the (x, y) coordinates
    BilinearInterpolator(
        const Vec<double>& xcoordinates,
        const Vec<double>& ycoordinates,
        const Fn<double(double, double)>& function);

    /// Set the x-coordinates of the interpolation.
    auto setCoordinatesX(const Vec<double>& xcoordinates) -> void;

    /// Set the y-coordinates of the interpolation.
    auto setCoordinatesY(const Vec<double>& ycoordinates) -> void;

    /// Set the data to be interpolated.
    auto setData(const Vec<double>& data) -> void;

    /// Set the data to be interpolated.
    auto setData(const Vec<Vec<double>>& data) -> void;

    /// Return the x-coordinates of the interpolation.
    auto xCoordinates() const -> const Vec<double>&;

    /// Return the y-coordinates of the interpolation.
    auto yCoordinates() const -> const Vec<double>&;

    /// Return the interpolation data.
    auto data() const -> const Vec<double>&;

    /// Check if the BilinearInterpolator instance is empty.
    auto empty() const -> bool;

    /// Calculate the interpolation at the provided (x, y) point.
    /// @param x The x-coordinate of the point
    /// @param y The y-coordinate of the point
    /// @return The interpolation of the data at (x, y) point
    auto operator()(const real& x, const real& y) const -> real;

private:
    /// The coordinates of the x and y points
    Vec<double> m_xcoordinates, m_ycoordinates;

    /// The interpolated data on every (x, y) point
    Vec<double> m_data;
};

/// Output a BilinearInterpolator instance
auto operator<<(std::ostream& out, const BilinearInterpolator& interpolator) -> std::ostream&;

} // namespace Reaktoro
