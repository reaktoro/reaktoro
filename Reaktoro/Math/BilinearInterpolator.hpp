// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

// C++ includes
#include <functional>
#include <iostream>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>

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
        const std::vector<real>& xcoordinates,
        const std::vector<real>& ycoordinates,
        const std::vector<real>& data);

    /// Construct a BilinearInterpolator instance with given function
    /// @param xcoordinates The x-coordinates for the interpolation
    /// @param ycoordinates The y-coordinates for the interpolation
    /// @param function The function to be interpolated over the (x, y) coordinates
    BilinearInterpolator(
        const std::vector<real>& xcoordinates,
        const std::vector<real>& ycoordinates,
        const std::function<real(real, real)>& function);

    /// setCoordinatesY the x-coordinates of the interpolation
    auto setCoordinatesX(const std::vector<real>& xcoordinates) -> void;

    /// Set the y-coordinates of the interpolation
    auto setCoordinatesY(const std::vector<real>& ycoordinates) -> void;

    /// Set the x-coordinates of the interpolation
    auto setData(const std::vector<real>& data) -> void;

    /// Return the x-coordinates of the interpolation
    auto xCoordinates() const -> const std::vector<real>&;

    /// Return the y-coordinates of the interpolation
    auto yCoordinates() const -> const std::vector<real>&;

    /// Return the interpolation data
    auto data() const -> const std::vector<real>&;

    /// Check if the BilinearInterpolator instance is empty
    auto empty() const -> bool;

    /// Calculate the interpolation at the provided (x, y) point
    /// @param x The x-coordinate of the point
    /// @param y The y-coordinate of the point
    /// @return The interpolation of the data at (x, y) point
    auto operator()(real x, real y) const -> real;

private:
    /// The coordinates of the x and y points
    std::vector<real> m_xcoordinates, m_ycoordinates;

    /// The interpolated data on every (x, y) point
    std::vector<real> m_data;
};

/// Output a BilinearInterpolator instance
auto operator<<(std::ostream& out, const BilinearInterpolator& interpolator) -> std::ostream&;

} // namespace Reaktoro
