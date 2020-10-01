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

#include "PhreeqcEngine.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Memoization.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcLegacy.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcThermo.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcWater.hpp>

namespace Reaktoro {
namespace PhreeqcUtils {

struct PhreeqcEngine::Impl
{
    /// The memoized function that calculates the thermodynamic and electrostatic properties of water just like in PHREEQC.
    Fn<PhreeqcWaterProps(real, real)> water_props_fn;

    /// Construct a default PhreeqcEngine::Impl object.
    Impl()
    {
        water_props_fn = PhreeqcUtils::waterProps;
        water_props_fn = memoizeLast(water_props_fn); // if given T,P same as last time, get cache result
    }

    /// Add missing standard molar volume property of a PHREEQC aqueous/exchange/surface species.
    auto addStandardVolume(StandardThermoProps& props, const PhreeqcSpecies* species, real T, real P) const -> void
    {
        const auto wprops = water_props_fn(T, P); // use memoized function in case of same T/P conditions as last time it was evaluated
        props.V0 = PhreeqcUtils::standardVolume(species, T, P, wprops); // in cm3/mol
        props.V0 *= cubicCentimeterToCubicMeter; // convert from cm3/mol to m3/mol
    }

    /// Add missing standard molar volume property of a PHREEQC gaseous/mineral species.
    auto addStandardVolume(StandardThermoProps& props, const PhreeqcPhase* phase, real T, real P) const -> void
    {
        props.V0 = PhreeqcUtils::standardVolume(phase, T, P); // in cm3/mol
        props.V0 *= cubicCentimeterToCubicMeter; // convert from cm3/mol to m3/mol
    }

    /// Add missing pressure correction to standard molar Gibbs energy and enthalpy of a PHREEQC aqueous/exchange/surface species.
    template<typename SpeciesType>
    auto addPressureCorrection(StandardThermoProps& props, const SpeciesType* species, real P) const -> void
    {
        const auto V0 = props.V0;
        const auto energy = PhreeqcUtils::pressureCorrectionEnergy(species, P, V0); // in J/mol
        props.G0 += energy;
        props.H0 += energy;
    }
};

PhreeqcEngine::PhreeqcEngine()
: pimpl(new Impl())
{}

auto PhreeqcEngine::addStandardVolume(StandardThermoProps& props, const PhreeqcSpecies* species, real T, real P) const -> void
{
    pimpl->addStandardVolume(props, species, T, P);
}

auto PhreeqcEngine::addStandardVolume(StandardThermoProps& props, const PhreeqcPhase* phase, real T, real P) const -> void
{
    pimpl->addStandardVolume(props, phase, T, P);
}

auto PhreeqcEngine::addPressureCorrection(StandardThermoProps& props, const PhreeqcSpecies* species, real P) const -> void
{
    pimpl->addPressureCorrection(props, species, P);
}

auto PhreeqcEngine::addPressureCorrection(StandardThermoProps& props, const PhreeqcPhase* phase, real P) const -> void
{
    pimpl->addPressureCorrection(props, phase, P);
}

} // namespace PhreeqcUtils
} // namespace Reaktoro
