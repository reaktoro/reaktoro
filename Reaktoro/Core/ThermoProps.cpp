// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "ThermoProps.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ThermoPropsPhase.hpp>

namespace Reaktoro {

ThermoProps::ThermoProps(const ChemicalSystem& system)
: sys(system)
{
    const auto numspecies = sys.species().size();
    const auto numphases = sys.phases().size();

    props.G0  = ArrayXr::Zero(numspecies);
    props.H0  = ArrayXr::Zero(numspecies);
    props.V0  = ArrayXr::Zero(numspecies);
    props.Cp0 = ArrayXr::Zero(numspecies);
    props.Cv0 = ArrayXr::Zero(numspecies);
}

ThermoProps::ThermoProps(const ChemicalSystem& system, const ThermoPropsData& data)
: sys(system), props(data)
{}

auto ThermoProps::system() const -> const ChemicalSystem&
{    return sys;
}

auto ThermoProps::data() const -> const ThermoPropsData&
{    return props;
}

auto ThermoProps::phase(Index idx) const -> ThermoPropsPhaseConstRef
{
    const auto phase = sys.phase(idx);
    const auto begin = sys.indexFirstSpeciesInPhase(idx);
    const auto size = phase.species().size();

    return ThermoPropsPhaseConstRef(phase, {
        props.T,
        props.P,
        props.G0.segment(begin, size),
        props.H0.segment(begin, size),
        props.V0.segment(begin, size),
        props.Cp0.segment(begin, size),
        props.Cv0.segment(begin, size),
    });
}

auto ThermoProps::phase(Index idx) -> ThermoPropsPhaseRef
{
    const auto phase = sys.phase(idx);
    const auto begin = sys.indexFirstSpeciesInPhase(idx);
    const auto size = phase.species().size();

    return ThermoPropsPhaseRef(phase, {
        props.T,
        props.P,
        props.G0.segment(begin, size),
        props.H0.segment(begin, size),
        props.V0.segment(begin, size),
        props.Cp0.segment(begin, size),
        props.Cv0.segment(begin, size),
    });
}

auto ThermoProps::temperature() const -> real
{
    return props.T;
}

auto ThermoProps::pressure() const -> real
{
    return props.P;
}

auto ThermoProps::standardGibbsEnergies() const -> ArrayXrConstRef
{
    return props.G0;
}

auto ThermoProps::standardEnthalpies() const -> ArrayXrConstRef
{
    return props.H0;
}

auto ThermoProps::standardVolumes() const -> ArrayXrConstRef
{
    return props.V0;
}

auto ThermoProps::standardEntropies() const -> ArrayXr
{
    return (props.H0 - props.G0)/props.T; // from G0 = H0 - T*S0
}

auto ThermoProps::standardInternalEnergies() const -> ArrayXr
{
    return props.H0 - props.P * props.V0; // from H0 = U0 + P*V0
}

auto ThermoProps::standardHelmholtzEnergies() const -> ArrayXr
{
    return props.G0 - props.P * props.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

auto ThermoProps::standardHeatCapacitiesConstP() const -> ArrayXrConstRef
{
    return props.Cp0;
}

auto ThermoProps::standardHeatCapacitiesConstV() const -> ArrayXrConstRef
{
    return props.Cv0;
}

} // namespace Reaktoro
