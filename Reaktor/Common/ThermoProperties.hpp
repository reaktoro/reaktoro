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

// Reaktor includes
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class  ThermoProperty;
struct ThermoPropertiesConstRow;
struct ThermoPropertiesRow;

/// Describe the thermodynamic properties and their partial temperature and pressure
/// derivatives of a collection of species or reactions.
/// ingroup Common
class ThermoProperties
{
public:
    /// Construct a ThermoProperties instance with given dimension
    /// @param nrows The number of rows of the vector quantities
    ThermoProperties(unsigned nrows);

    /// Construct a ThermoProperties instance
    /// @param val The values of the thermodynamic properties
    /// @param ddt The partial temperature derivatives of the thermodynamic properties
    /// @param ddp The partial pressure derivatives of the thermodynamic properties
    ThermoProperties(const Vector& val, const Vector& ddt, const Vector& ddp);

    /// Get the values of the thermodynamic properties
    auto val() const -> const Vector&;

    /// Get the partial temperature derivatives of the thermodynamic properties
    auto ddt() const -> const Vector&;

    /// Get the partial pressure derivatives of the thermodynamic properties
    auto ddp() const -> const Vector&;

    /// Get a reference of a row of this ThermoProperties instance
    auto row(unsigned irow) -> ThermoPropertiesRow;

    /// Get a const reference of a row of this ThermoProperties instance
    auto row(unsigned irow) const -> ThermoPropertiesConstRow;

private:
    /// The values of the thermodynamic properties
    Vector m_val;

    /// The partial temperature derivatives of the thermodynamic properties
    Vector m_ddt;

    /// The partial pressure derivatives of the thermodynamic properties
    Vector m_ddp;
};

/// An auxiliary type for the representation of the view of a row of a ThermoProperties instance
struct ThermoPropertiesRow
{
    ThermoPropertiesRow(const ThermoProperties& vector, unsigned irow);
    auto operator=(const ThermoProperty& property) -> ThermoPropertiesRow&;
    VectorRow val;
    VectorRow ddt;
    VectorRow ddp;
};

/// An auxiliary type for the representation of the const view of a row of a ThermoProperties instance
struct ThermoPropertiesConstRow
{
    ThermoPropertiesConstRow(const ThermoProperties& properties, unsigned irow);
    const VectorRow val;
    const VectorRow ddt;
    const VectorRow ddp;
};

/// Compares two ThermoProperties instances for equality
auto operator==(const ThermoProperties& l, const ThermoProperties& r) -> bool;

} // namespace Reaktor
