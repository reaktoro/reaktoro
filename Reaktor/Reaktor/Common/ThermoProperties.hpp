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

/// Describe the thermodynamic properties and their partial temperature and pressure
/// derivatives of a collection of species or reactions.
/// ingroup Common
struct ThermoProperties
{
    /// Construct a ThermoProperties instance
    /// @param val The values of the thermodynamic properties of the species
    /// @param ddt The partial temperature derivatives of the thermodynamic properties of the species
    /// @param ddp The partial pressure derivatives of the thermodynamic properties of the species
    ThermoProperties(const Vector& val, const Vector& ddt, const Vector& ddp);

    /// Construct a ThermoProperties instance from rvalues
    /// @param val The values of the thermodynamic properties of the species
    /// @param ddt The partial temperature derivatives of the thermodynamic properties of the species
    /// @param ddp The partial pressure derivatives of the thermodynamic properties of the species
    ThermoProperties(Vector&& val, Vector&& ddt, Vector&& ddp);

    /// The values of the thermodynamic property of every species
    const Vector val;

    /// The partial derivatives of the thermodynamic property of every species w.r.t. temperature (in units of K)
    const Vector ddt;

    /// The partial derivatives of the thermodynamic property of every species w.r.t. pressure (in units of Pa)
    const Vector ddp;
};

}  // namespace Reaktor
