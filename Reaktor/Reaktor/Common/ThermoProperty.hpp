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

namespace Reaktor {

/// Describe a thermodynamic property value and its partial derivatives w.r.t. temperature and pressure
struct ThermoProperty
{
    /// The value of the thermodynamic property
    double val = 0.0;

    /// The partial derivative of the thermodynamic property w.r.t. temperature (in units of K)
    double ddt = 0.0;

    /// The partial derivative of the thermodynamic property w.r.t. pressure (in units of Pa)
    double ddp = 0.0;
};

/// Describe the function signature of a thermodynamic property function
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
/// @return A ThermoProperty instance with the thermodynamic property of a species
/// @ see ThermoProperty
typedef std::function<ThermoProperty(double T, double P)> ThermoPropertyFunction;

}  // namespace Reaktor
