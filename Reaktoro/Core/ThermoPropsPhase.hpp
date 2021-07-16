// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/ArrayStream.hpp>
#include <Reaktoro/Common/AutoDiff.hpp>
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

/// The base type for primary standard thermodynamic property data of a phase from which others are computed.
template<typename Real, typename Array>
struct ThermoPropsPhaseBaseData
{
    /// The temperature of the phase (in K).
    Real T;

    /// The pressure of the phase (in Pa).
    Real P;

    /// The standard molar Gibbs energies of formation of the species in the phase (in J/mol)
    Array G0;

    /// The standard molar enthalpies of formation of the species in the phase (in J/mol)
    Array H0;

    /// The standard molar volumes of the species in the phase (in m3/mol)
    Array V0;

    /// The standard molar isobaric heat capacities of the species in the phase (in J/(mol·K))
    Array Cp0;

    /// The standard molar isochoric heat capacities of the species in the phase (in J/(mol·K))
    Array Cv0;

    /// Assign a ThermoPropsPhaseBaseData object to this.
    template<typename RX, typename AX>
    auto operator=(const ThermoPropsPhaseBaseData<RX, AX>& other)
    {
        T   = other.T;
        P   = other.P;
        G0  = other.G0;
        H0  = other.H0;
        V0  = other.V0;
        Cp0 = other.Cp0;
        Cv0 = other.Cv0;
        return *this;
    }

    /// Convert this ThermoPropsPhaseBaseData object into another.
    template<typename RX, typename AX>
    operator ThermoPropsPhaseBaseData<RX, AX>()
    {
        return { T, P, G0, H0, V0, Cp0, Cv0 };
    }

    /// Convert this ThermoPropsPhaseBaseData object into another.
    template<typename RX, typename AX>
    operator ThermoPropsPhaseBaseData<RX, AX>() const
    {
        return { T, P, G0, H0, V0, Cp0, Cv0 };
    }

    /// Assign the given array data to this ThermoPropsPhaseBaseData object.
    auto operator=(const ArrayStream<Real>& array)
    {
        array.to(T, P, G0, H0, V0, Cp0, Cv0);
        return *this;
    }

    /// Convert this ThermoPropsPhaseBaseData object into an array.
    operator ArrayStream<Real>() const
    {
        return {T, P, G0, H0, V0, Cp0, Cv0};
    }
};

/// The primary standard thermodynamic property data of a phase from which others are computed.
using ThermoPropsPhaseData = ThermoPropsPhaseBaseData<real, ArrayXr>;

/// The primary standard thermodynamic property data of a phase from which others are computed.
using ThermoPropsPhaseDataRef = ThermoPropsPhaseBaseData<real&, ArrayXrRef>;

/// The primary standard thermodynamic property data of a phase from which others are computed.
using ThermoPropsPhaseDataConstRef = ThermoPropsPhaseBaseData<const real&, ArrayXrConstRef>;

/// The type of functions that computes the primary standard thermodynamic property data of a phase.
using ThermoPropsPhaseFn = Fn<void(ThermoPropsPhaseDataRef, const real&, const real&, ArrayXrConstRef)>;

/// The base type for standard thermodynamic properties of a phase and its species.
template<typename Real, typename Array>
class ThermoPropsPhaseBase
{
public:
    /// Construct a ThermoPropsPhaseBase instance.
    explicit ThermoPropsPhaseBase(const Phase& phase)
    : mphase(phase)
    {
        const auto numspecies = phase.species().size();

        mdata.G0   = ArrayXr::Zero(numspecies);
        mdata.H0   = ArrayXr::Zero(numspecies);
        mdata.V0   = ArrayXr::Zero(numspecies);
        mdata.Cp0  = ArrayXr::Zero(numspecies);
        mdata.Cv0  = ArrayXr::Zero(numspecies);
    }

    /// Construct a ThermoPropsPhaseBase instance.
    ThermoPropsPhaseBase(const Phase& phase, const ThermoPropsPhaseBaseData<Real, Array>& data)
    : mphase(phase), mdata(data)
    {}

    /// Construct a ThermoPropsPhaseBase instance.
    template<typename RX, typename AX>
    ThermoPropsPhaseBase(ThermoPropsPhaseBase<RX, AX>& other)
    : mphase(other.mphase), mdata(other.mdata)
    {}

    /// Construct a ThermoPropsPhaseBase instance.
    template<typename RX, typename AX>
    ThermoPropsPhaseBase(const ThermoPropsPhaseBase<RX, AX>& other)
    : mphase(other.mphase), mdata(other.mdata)
    {}

    /// Update the standard thermodynamic properties of the phase.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    auto update(const real& T, const real& P)
    {
        // Check if this update call can be skipped if T, P conditions remain the same
        if(T == mdata.T && P == mdata.P)
            return;

        mdata.T = T;
        mdata.P = P;

        auto& G0   = mdata.G0;
        auto& H0   = mdata.H0;
        auto& V0   = mdata.V0;
        auto& Cp0  = mdata.Cp0;
        auto& Cv0  = mdata.Cv0;

        const auto& species = mphase.species();
        const auto size = species.size();

        assert(G0.size()   == size);
        assert(H0.size()   == size);
        assert(V0.size()   == size);
        assert(Cp0.size()  == size);
        assert(Cv0.size()  == size);

        // Compute the standard thermodynamic properties of the species in the phase.
        StandardThermoProps aux;
        for(auto i = 0; i < size; ++i)
        {
            const auto& standard_thermo_model = mphase.species(i).standardThermoModel();
            aux = standard_thermo_model ? standard_thermo_model(T, P) : StandardThermoProps{};
            G0[i]  = aux.G0;
            H0[i]  = aux.H0;
            V0[i]  = aux.V0;
            Cp0[i] = aux.Cp0;
            Cv0[i] = aux.Cv0;
        }
    }

    /// Update the standard thermodynamic properties of the phase.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param wrtvar The variable with respect to automatic differentiation should be carried out.
    auto update(const real& T, const real& P, Wrt<real&> wrtvar)
    {
        autodiff::seed(wrtvar);
        update(T, P);
        autodiff::unseed(wrtvar);
    }

    /// Return the underlying Phase object.
    auto phase() const -> const Phase&
    {
        return mphase;
    }

    /// Return the primary standard thermodynamic property data of the phase from which others are calculated.
    auto data() const -> const ThermoPropsPhaseBaseData<Real, Array>&
    {
        return mdata;
    }

    /// Return the temperature of the phase (in K).
    auto temperature() const -> real
    {
        return mdata.T;
    }

    /// Return the pressure of the phase (in Pa).
    auto pressure() const -> real
    {
        return mdata.P;
    }

    /// Return the standard partial molar volumes of the species (in m3/mol).
    auto standardVolumes() const -> ArrayXrConstRef
    {
        return mdata.V0;
    }

    /// Return the standard partial molar Gibbs energies of formation of the species (in J/mol).
    auto standardGibbsEnergies() const -> ArrayXrConstRef
    {
        return mdata.G0;
    }

    /// Return the standard partial molar enthalpies of formation of the species (in J/mol).
    auto standardEnthalpies() const -> ArrayXrConstRef
    {
        return mdata.H0;
    }

    /// Return the standard partial molar entropies of formation of the species (in J/(mol*K)).
    auto standardEntropies() const -> ArrayXr
    {
        return (mdata.H0 - mdata.G0)/mdata.T; // from G0 = H0 - T*S0
    }

    /// Return the standard partial molar internal energies of formation of the species (in J/mol).
    auto standardInternalEnergies() const -> ArrayXr
    {
        return mdata.H0 - mdata.P * mdata.V0; // from H0 = U0 + P*V0
    }

    /// Return the standard partial molar Helmholtz energies of formation of the species (in J/mol).
    auto standardHelmholtzEnergies() const -> ArrayXr
    {
        return mdata.G0 - mdata.P * mdata.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
    }

    /// Return the standard partial molar isobaric heat capacities of the species (in J/(mol*K)).
    auto standardHeatCapacitiesConstP() const -> ArrayXrConstRef
    {
        return mdata.Cp0;
    }

    /// Return the standard partial molar isochoric heat capacities of the species (in J/(mol*K)).
    auto standardHeatCapacitiesConstV() const -> ArrayXrConstRef
    {
        return mdata.Cv0;
    }

    /// Assign the given array data to this ThermoPropsPhaseBase object.
    auto operator=(const ArrayStream<Real>& array)
    {
        mdata = array;
        return *this;
    }

    /// Convert this ThermoPropsPhaseBase object into an array.
    operator ArrayStream<Real>() const
    {
        return mdata;
    }

    // Ensure other ThermoPropsPhaseBase types are friend among themselves.
    template<typename RX, typename AX>
    friend class ThermoPropsPhaseBase;

private:
    /// The phase associated with these primary standard thermodynamic properties.
    Phase mphase;

    /// The primary standard thermodynamic property data of the phase from which others are calculated.
    ThermoPropsPhaseBaseData<Real, Array> mdata;
};

/// The standard thermodynamic properties of a phase and its species.
using ThermoPropsPhase = ThermoPropsPhaseBase<real, ArrayXr>;

/// The non-const view to the standard thermodynamic properties of a phase and its species.
using ThermoPropsPhaseRef = ThermoPropsPhaseBase<real&, ArrayXrRef>;

/// The const view to the standard thermodynamic properties of a phase and its species.
using ThermoPropsPhaseConstRef = ThermoPropsPhaseBase<const real&, ArrayXrConstRef>;

} // namespace Reaktoro
