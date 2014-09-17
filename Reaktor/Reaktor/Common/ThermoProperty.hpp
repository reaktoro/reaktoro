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

namespace Reaktor {

/// Describe a thermodynamic property value and its partial derivatives w.r.t. temperature and pressure
class ThermoProperty
{
public:
    /// Construct a default ThermoProperty instance
    ThermoProperty();

    /// Construct a ThermoProperty instance
    /// @param val The value of the thermodynamic property
    /// @param ddt The partial temperature derivative of the thermodynamic property
    /// @param ddp The partial pressure derivative of the thermodynamic property
    ThermoProperty(double val, double ddt, double ddp);

    /// Get the value of the thermodynamic property
    auto val() const -> double;

    /// Get the partial temperature derivative of the thermodynamic property
    auto ddt() const -> double;

    /// Get the partial pressure derivative of the thermodynamic property
    auto ddp() const -> double;

private:
    /// The value of the thermodynamic property
    double m_val = 0.0;

    /// The partial temperature derivative of the thermodynamic property
    double m_ddt = 0.0;

    /// The partial pressure derivative of the thermodynamic property
    double m_ddp = 0.0;
};

/// Describe the function signature of a thermodynamic property function
/// @param T The temperature (in units of K)
/// @param P The pressure (in units of Pa)
/// @return A ThermoProperty instance with the thermodynamic property of a species
/// @ see ThermoProperty
typedef std::function<ThermoProperty(double T, double P)> ThermoPropertyFunction;

}  // namespace Reaktor
