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
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

/// The base type for primary standard thermodynamic property data of a phase.
template<typename Real, typename Array>
struct PhaseThermoPropsBaseData;

/// The primary standard thermodynamic property data of a phase.
using PhaseThermoPropsData = PhaseThermoPropsBaseData<real, ArrayXr>;

/// The primary standard thermodynamic property data of a phase.
using PhaseThermoPropsDataRef = PhaseThermoPropsBaseData<real&, ArrayXrRef>;

/// The primary standard thermodynamic property data of a phase.
using PhaseThermoPropsDataConstRef = PhaseThermoPropsBaseData<const real&, ArrayXrConstRef>;

/// The base type for standard thermodynamic properties of a phase and its species.
template<typename R, typename A>
class PhaseThermoPropsBase;

/// The standard thermodynamic properties of a phase and its species.
using PhaseThermoProps = PhaseThermoPropsBase<real, ArrayXr>;

/// The non-const view to the standard thermodynamic properties of a phase and its species.
using PhaseThermoPropsRef = PhaseThermoPropsBase<real&, ArrayXrRef>;

/// The const view to the standard thermodynamic properties of a phase and its species.
using PhaseThermoPropsConstRef = PhaseThermoPropsBase<const real&, ArrayXrConstRef>;

/// The base type for primary standard thermodynamic property data of a phase.
template<typename Real, typename Array>
struct PhaseThermoPropsBaseData
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
class PhaseThermoPropsBase
{
public:
    /// Construct a PhaseThermoPropsBase instance.
    explicit PhaseThermoPropsBase(const Phase& phase);

    /// Construct a PhaseThermoPropsBase instance.
    PhaseThermoPropsBase(const Phase& phase, const PhaseThermoPropsBaseData<R, A>& data);

    /// Construct a PhaseThermoPropsBase instance.
    template<typename RX, typename AX>
    PhaseThermoPropsBase(PhaseThermoPropsBase<RX, AX>& props);

    /// Construct a PhaseThermoPropsBase instance.
    template<typename RX, typename AX>
    PhaseThermoPropsBase(const PhaseThermoPropsBase<RX, AX>& props);

    /// Update the standard thermodynamic properties of the phase and its species.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    auto update(real T, real P);

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

    // Ensure other PhaseThermoPropsBase types are friend among themselves.
    template<typename RX, typename AX>
    friend class PhaseThermoPropsBase;

private:
    /// The phase associated with these primary standard thermodynamic properties.
    Phase phase;

    /// The primary standard thermodynamic property data of the phase from which others are calculated.
    PhaseThermoPropsBaseData<R, A> props;
};

template<typename Real, typename Array>
PhaseThermoPropsBase<Real, Array>::PhaseThermoPropsBase(const Phase& phase)
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
auto PhaseThermoPropsBase<Real, Array>::update(real Tnew, real Pnew)
{
    auto& T   = props.T;
    auto& P   = props.P;
    auto& G0  = props.G0;
    auto& H0  = props.H0;
    auto& V0  = props.V0;
    auto& Cp0 = props.Cp0;
    auto& Cv0 = props.Cv0;

    const auto& species = phase.species();
    const auto size = species.size();

    assert(G0.size()   == size);
    assert(H0.size()   == size);
    assert(V0.size()   == size);
    assert(Cp0.size()  == size);
    assert(Cv0.size()  == size);

    const auto updatingT = T != Tnew;
    const auto updatingP = P != Pnew;

    // The function that evaluates the standard thermodynamic properties of the species in the phase.
    auto evalStandardThermoProps = [&]()
    {
        StandardThermoProps aux;
        for(auto i = 0; i < size; ++i)
        {
            assert(phase.species(i).standardThermoPropsFn());
            aux = phase.species(i).standardThermoPropsFn()(T, P);
            G0[i]  = aux.G0;
            H0[i]  = aux.H0;
            V0[i]  = aux.V0;
            Cp0[i] = aux.Cp0;
            Cv0[i] = aux.Cv0;
        }
    };

    if(updatingT || updatingP)
    {
        T = Tnew;
        P = Pnew;
        evalStandardThermoProps();
    }
}

template<typename Real, typename Array>
PhaseThermoPropsBase<Real, Array>::PhaseThermoPropsBase(const Phase& phase, const PhaseThermoPropsBaseData<Real, Array>& data)
: phase(phase), props(data)
{}

template<typename Real, typename Array>
template<typename RX, typename AX>
PhaseThermoPropsBase<Real, Array>::PhaseThermoPropsBase(PhaseThermoPropsBase<RX, AX>& other)
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
PhaseThermoPropsBase<Real, Array>::PhaseThermoPropsBase(const PhaseThermoPropsBase<RX, AX>& other)
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
auto PhaseThermoPropsBase<Real, Array>::data() const
{
    return props;
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::temperature() const
{
    return props.T;
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::pressure() const
{
    return props.P;
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::standardGibbsEnergies() const
{
    return props.G0;
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::standardEnthalpies() const
{
    return props.H0;
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::standardVolumes() const
{
    return props.V0;
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::standardEntropies() const
{
    return (props.H0 - props.G0)/props.T; // from G0 = H0 - T*S0
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::standardInternalEnergies() const
{
    return props.H0 - props.P * props.V0; // from H0 = U0 + P*V0
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::standardHelmholtzEnergies() const
{
    return props.G0 - props.P * props.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::standardHeatCapacitiesConstP() const
{
    return props.Cp0;
}

template<typename Real, typename Array>
auto PhaseThermoPropsBase<Real, Array>::standardHeatCapacitiesConstV() const
{
    return props.Cv0;
}

} // namespace Reaktoro
