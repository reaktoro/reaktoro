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
#include <functional>
#include <iostream>
#include <vector>

namespace Reaktor {

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
        const std::vector<double>& xcoordinates,
        const std::vector<double>& ycoordinates,
        const std::vector<double>& data);

    /// Construct a BilinearInterpolator instance with given function
    /// @param xcoordinates The x-coordinates for the interpolation
    /// @param ycoordinates The y-coordinates for the interpolation
    /// @param function The function to be interpolated over the (x, y) coordinates
    BilinearInterpolator(
        const std::vector<double>& xcoordinates,
        const std::vector<double>& ycoordinates,
        const std::function<double(double, double)>& function);

    /// setCoordinatesY the x-coordinates of the interpolation
    auto setCoordinatesX(const std::vector<double>& xcoordinates) -> void;

    /// Set the y-coordinates of the interpolation
    auto setCoordinatesY(const std::vector<double>& ycoordinates) -> void;

    /// Set the x-coordinates of the interpolation
    auto setData(const std::vector<double>& data) -> void;

    /// Get the x-coordinates of the interpolation
    auto xCoodinates() const -> const std::vector<double>&;

    /// Get the y-coordinates of the interpolation
    auto yCoodinates() const -> const std::vector<double>&;

    /// Get the interpolation data
    auto data() const -> const std::vector<double>&;

    /// Check if the BilinearInterpolator instance was initialised
    auto initialised() const -> bool;

    /// Calculate the interpolation at the provided (x, y) point
    /// @param x The x-coordinate of the point
    /// @param y The y-coordinate of the point
    /// @return The interpolation of the data at (x, y) point
    auto operator()(double x, double y) const -> double;

private:
    /// The coordinates of the x and y points
    std::vector<double> xcoordinates$, ycoordinates$;

    /// The interpolated data on every (x, y) point
    std::vector<double> data$;
};

/// Output a BilinearInterpolator instance
auto operator<<(std::ostream& out, const BilinearInterpolator& interpolator) -> std::ostream&;

} /* namespace Reaktor */
