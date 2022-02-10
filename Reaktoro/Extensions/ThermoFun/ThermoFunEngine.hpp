// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/StandardThermoProps.hpp>

// Forward declaration of ThermoFun classes
namespace ThermoFun { class Database; class Substance; }

namespace Reaktoro {

/// The class used for standard thermodynamic property calculations based on ThermoFun.
/// @ingroup ThermoFunExtension
class ThermoFunEngine
{
public:
    /// Construct a ThermoFunEngine object with given database.
    ThermoFunEngine(const ThermoFun::Database& database);

    /// Return the ThermoFun::Database object.
    auto database() const -> const ThermoFun::Database&;

    /// Return the standard thermodynamic properties of a chemical species with given name.
    auto props(const real& T, const real& P, const String& species) const -> StandardThermoProps;

    /// Return the standard thermodynamic properties of a chemical species with given ThermoFun::Substance object.
    auto props(const real& T, const real& P, const ThermoFun::Substance& substance) const -> StandardThermoProps;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

} // namespace Reaktoro
