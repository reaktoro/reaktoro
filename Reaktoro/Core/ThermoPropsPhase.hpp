// Reaktoro is a unified framework for modeling chemically reactive phases.
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
#include <Reaktoro/Core/ThermoPropsPhase.fwd.hpp>
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

/// The base type for primary standard thermodynamic property data of a phase.
template<typename Real, typename Array>
struct ThermoPropsPhaseBaseData
{
    /// The temperature of the phase (in K).
    Real T;

    /// The pressure of the phase (in Pa).
    Real P;

    /// The standard molar Gibbs energies of the species in the phase (in J/mol)
    Array G0;

    /// The standard molar enthalpies of the species in the phase (in J/mol)
    Array H0;

    /// The standard molar volumes of the species in the phase (in m3/mol)
    Array V0;

    /// The standard molar isobaric heat capacities of the species in the phase (in J/(mol·K))
    Array Cp0;

    /// The standard molar isochoric heat capacities of the species in the phase (in J/(mol·K))
    Array Cv0;
};

/// The base type for standard thermodynamic properties of a phase and its species.
template<typename R, typename A>
class ThermoPropsPhaseBase
{
public:
    /// Construct a ThermoPropsPhaseBase instance.
    explicit ThermoPropsPhaseBase(const Phase& phase);

    /// Construct a ThermoPropsPhaseBase instance.
    ThermoPropsPhaseBase(const Phase& phase, const ThermoPropsPhaseBaseData<R, A>& data);

    /// Construct a ThermoPropsPhaseBase instance.
    template<typename RX, typename AX>
    ThermoPropsPhaseBase(ThermoPropsPhaseBase<RX, AX>& props);

    /// Construct a ThermoPropsPhaseBase instance.
    template<typename RX, typename AX>
    ThermoPropsPhaseBase(const ThermoPropsPhaseBase<RX, AX>& props);

    /// Return the primary chemical property data of the phase from which others are calculated.
    auto data() const;

    /// Return the temperature of the phase (in K).
    auto temperature() const;

    /// Return the pressure of the phase (in Pa).
    auto pressure() const;

    /// Return the standard partial molar Gibbs energies of the species (in J/mol).
    auto standardGibbsEnergies() const;

    /// Return the standard partial molar enthalpies of the species (in J/mol).
    auto standardEnthalpies() const;

    /// Return the standard partial molar volumes of the species (in m3/mol).
    auto standardVolumes() const;

    /// Return the standard partial molar entropies of the species (in J/(mol*K)).
    auto standardEntropies() const;

    /// Return the standard partial molar internal energies of the species (in J/mol).
    auto standardInternalEnergies() const;

    /// Return the standard partial molar Helmholtz energies of the species (in J/mol).
    auto standardHelmholtzEnergies() const;

    /// Return the standard partial molar isobaric heat capacities of the species (in J/(mol*K)).
    auto standardHeatCapacitiesConstP() const;

    /// Return the standard partial molar isochoric heat capacities of the species (in J/(mol*K)).
    auto standardHeatCapacitiesConstV() const;

    // Ensure other ThermoPropsPhaseBase types are friend among themselves.
    template<typename RX, typename AX>
    friend class ThermoPropsPhaseBase;

private:
    /// The phase associated with these primary standard thermodynamic properties.
    Phase phase;

    /// The primary standard thermodynamic property data of the phase from which others are calculated.
    ThermoPropsPhaseBaseData<R, A> props;
};

template<typename Real, typename Array>
ThermoPropsPhaseBase<Real, Array>::ThermoPropsPhaseBase(const Phase& phase)
: phase(phase)
{
    const auto numspecies = phase.species().size();

    props.G0   = ArrayXr::Zero(numspecies);
    props.H0   = ArrayXr::Zero(numspecies);
    props.V0   = ArrayXr::Zero(numspecies);
    props.Cp0  = ArrayXr::Zero(numspecies);
    props.Cv0  = ArrayXr::Zero(numspecies);
}

template<typename Real, typename Array>
ThermoPropsPhaseBase<Real, Array>::ThermoPropsPhaseBase(const Phase& phase, const ThermoPropsPhaseBaseData<Real, Array>& data)
: phase(phase), props(data)
{}

template<typename Real, typename Array>
template<typename RX, typename AX>
ThermoPropsPhaseBase<Real, Array>::ThermoPropsPhaseBase(ThermoPropsPhaseBase<RX, AX>& other)
: phase(other.phase), props{
    other.props.T,
    other.props.P,
    other.props.G0,
    other.props.H0,
    other.props.V0,
    other.props.Cp0,
    other.props.Cv0 }
{}

template<typename Real, typename Array>
template<typename RX, typename AX>
ThermoPropsPhaseBase<Real, Array>::ThermoPropsPhaseBase(const ThermoPropsPhaseBase<RX, AX>& other)
: phase(other.phase), props{
    other.props.T,
    other.props.P,
    other.props.G0,
    other.props.H0,
    other.props.V0,
    other.props.Cp0,
    other.props.Cv0 }
{}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::data() const
{
    return props;
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::temperature() const
{
    return props.T;
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::pressure() const
{
    return props.P;
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::standardGibbsEnergies() const
{
    return props.G0;
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::standardEnthalpies() const
{
    return props.H0;
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::standardVolumes() const
{
    return props.V0;
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::standardEntropies() const
{
    return (props.H0 - props.G0)/props.T; // from G0 = H0 - T*S0
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::standardInternalEnergies() const
{
    return props.H0 - props.P * props.V0; // from H0 = U0 + P*V0
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::standardHelmholtzEnergies() const
{
    return props.G0 - props.P * props.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::standardHeatCapacitiesConstP() const
{
    return props.Cp0;
}

template<typename Real, typename Array>
auto ThermoPropsPhaseBase<Real, Array>::standardHeatCapacitiesConstV() const
{
    return props.Cv0;
}

} // namespace Reaktoro
