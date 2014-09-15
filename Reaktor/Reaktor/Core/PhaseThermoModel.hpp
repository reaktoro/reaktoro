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

// Reaktor includes
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>

namespace Reaktor {

/// Describe the thermodynamic model of a phase
struct PhaseThermoModel
{
    /// The activity function of the phase
    std::function<ThermoVector(double, double, const Vector&)> activity;

    /// The concentration function of the phase
    std::function<ThermoVector(const Vector&)> concentration;

    /// The density function of the phase (in units of kg/m3)
    std::function<ThermoScalar(double, double, const Vector&)> density;
};

}  // namespace Reaktor


