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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/StandardThermoProps.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp> // ***BECAUSE PhreeqcLegacy.hpp IS INCLUDED HERE, MAKE SURE THIS HEADER FILE IS NOT EXPORTED!***

namespace Reaktoro {
namespace PhreeqcUtils {

/// The class used for standard thermodynamic property calculations based on SUPCRT.
/// @ingroup PhreeqcExtension
class PhreeqcEngine
{
public:
    /// Construct a default PhreeqcEngine object.
    PhreeqcEngine();

    /// Add missing standard molar volume property of a PHREEQC aqueous/exchange/surface species.
    auto addStandardVolume(StandardThermoProps& props, const PhreeqcSpecies* species, real T, real P) const -> void;

    /// Add missing standard molar volume property of a PHREEQC gaseous/mineral species.
    auto addStandardVolume(StandardThermoProps& props, const PhreeqcPhase* phase, real T, real P) const -> void;

    /// Add missing pressure correction to standard molar Gibbs energy and enthalpy of a PHREEQC aqueous/exchange/surface species.
    auto addPressureCorrection(StandardThermoProps& props, const PhreeqcSpecies* species, real P) const -> void;

    /// Add missing pressure correction to standard molar Gibbs energy and enthalpy of a PHREEQC gaseous/mineral species.
    auto addPressureCorrection(StandardThermoProps& props, const PhreeqcPhase* phase, real P) const -> void;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

} // namespace PhreeqcUtils
} // namespace Reaktoro
