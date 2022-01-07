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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

/// The base type for primary chemical property data of a phase from which others are computed.
template<typename Real, typename Array>
struct ChemicalPropsPhaseBaseData
{
    /// The temperature of the phase (in K).
    Real T;

    /// The pressure of the phase (in Pa).
    Real P;

    /// The amounts of each species in the phase (in mol).
    Array n;

    /// The amount of the phase as the sum of species amounts (in mol).
    Real nsum;

    /// The mass of the phase as the sum of species masses (in kg).
    Real msum;

    /// The mole fractions of the species in the phase (in mol/mol).
    Array x;

    /// The standard molar Gibbs energies of formation of the species in the phase (in J/mol)
    Array G0;

    /// The standard molar enthalpies of formation of the species in the phase (in J/mol)
    Array H0;

    /// The standard molar volumes of the species in the phase (in m³/mol)
    Array V0;

    /// The temperature derivative of the standard molar volumes of the species in the phase (in m³/(mol·K)).
    Array VT0;

    /// The pressure derivative of the standard molar volumes of the species in the phase (in m³/(mol·Pa)).
    Array VP0;

    /// The standard molar isobaric heat capacities of the species in the phase (in J/(mol·K))
    Array Cp0;

    /// The corrective molar volume of the phase (in m³/mol).
    Real Vx;

    /// The temperature derivative of the corrective molar volume at constant pressure (in m³/(mol·K)).
    Real VxT;

    /// The pressure derivative of the corrective molar volume at constant temperature (in m³/(mol·Pa)).
    Real VxP;

    /// The corrective molar Gibbs energy of the phase (in units of J/mol).
    Real Gx;

    /// The corrective molar enthalpy of the phase (in units of J/mol).
    Real Hx;

    /// The corrective molar isobaric heat capacity of the phase (in units of J/(mol·K)).
    Real Cpx;

    /// The activity coefficients (natural log) of the species in the phase.
    Array ln_g;

    /// The activities (natural log) of the species in the phase.
    Array ln_a;

    /// The chemical potentials of the species in the phase.
    Array u;

    /// Assign a ChemicalPropsPhaseBaseData object to this.
    template<typename RX, typename AX>
    auto operator=(const ChemicalPropsPhaseBaseData<RX, AX>& other)
    {
        T    = other.T;
        P    = other.P;
        n    = other.n;
        nsum = other.nsum;
        msum = other.msum;
        x    = other.x;
        G0   = other.G0;
        H0   = other.H0;
        V0   = other.V0;
        VT0  = other.VT0;
        VP0  = other.VP0;
        Cp0  = other.Cp0;
        Vx  = other.Vx;
        VxT = other.VxT;
        VxP = other.VxP;
        Gx  = other.Gx;
        Hx  = other.Hx;
        Cpx = other.Cpx;
        ln_g = other.ln_g;
        ln_a = other.ln_a;
        u    = other.u;
        return *this;
    }

    /// Convert this ChemicalPropsPhaseBaseData object into another.
    template<typename RX, typename AX>
    operator ChemicalPropsPhaseBaseData<RX, AX>()
    {
        return { T, P, n, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Gx, Hx, Cpx, ln_g, ln_a, u };
    }

    /// Convert this ChemicalPropsPhaseBaseData object into another.
    template<typename RX, typename AX>
    operator ChemicalPropsPhaseBaseData<RX, AX>() const
    {
        return { T, P, n, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Gx, Hx, Cpx, ln_g, ln_a, u };
    }

    /// Assign the given array data to this ChemicalPropsPhaseBaseData object.
    auto operator=(const ArrayStream<Real>& array)
    {
        array.to(T, P, n, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Gx, Hx, Cpx, ln_g, ln_a, u);
        return *this;
    }

    /// Convert this ChemicalPropsPhaseBaseData object into an array.
    operator ArrayStream<Real>() const
    {
        return { T, P, n, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Gx, Hx, Cpx, ln_g, ln_a, u };
    }
};

/// The primary chemical property data of a phase from which others are computed.
using ChemicalPropsPhaseData = ChemicalPropsPhaseBaseData<real, ArrayXr>;

/// The primary chemical property data of a phase from which others are computed.
using ChemicalPropsPhaseDataRef = ChemicalPropsPhaseBaseData<real&, ArrayXrRef>;

/// The primary chemical property data of a phase from which others are computed.
using ChemicalPropsPhaseDataConstRef = ChemicalPropsPhaseBaseData<const real&, ArrayXrConstRef>;

/// The type of functions that computes the primary chemical property data of a phase.
using ChemicalPropsPhaseFn = Fn<void(ChemicalPropsPhaseDataRef, const real&, const real&, ArrayXrConstRef)>;

/// The base type for chemical properties of a phase and its species.
template<typename Real, typename Array>
class ChemicalPropsPhaseBase
{
public:
    /// Construct a ChemicalPropsPhaseBase instance.
    explicit ChemicalPropsPhaseBase(const Phase& phase)
    : mphase(phase)
    {
        const auto N = phase.species().size();

        mdata.n    = ArrayXr::Zero(N);
        mdata.x    = ArrayXr::Zero(N);
        mdata.G0   = ArrayXr::Zero(N);
        mdata.H0   = ArrayXr::Zero(N);
        mdata.V0   = ArrayXr::Zero(N);
        mdata.VT0  = ArrayXr::Zero(N);
        mdata.VP0  = ArrayXr::Zero(N);
        mdata.Cp0  = ArrayXr::Zero(N);
        mdata.ln_g = ArrayXr::Zero(N);
        mdata.ln_a = ArrayXr::Zero(N);
        mdata.u    = ArrayXr::Zero(N);
    }

    /// Construct a ChemicalPropsPhaseBase instance.
    ChemicalPropsPhaseBase(const Phase& phase, const ChemicalPropsPhaseBaseData<Real, Array>& data)
    : mphase(phase), mdata(data)
    {}

    /// Construct a ChemicalPropsPhaseBase instance.
    template<typename RX, typename AX>
    ChemicalPropsPhaseBase(ChemicalPropsPhaseBase<RX, AX>& other)
    : mphase(other.mphase), mdata(other.mdata)
    {}

    /// Construct a ChemicalPropsPhaseBase instance.
    template<typename RX, typename AX>
    ChemicalPropsPhaseBase(const ChemicalPropsPhaseBase<RX, AX>& other)
    : mphase(other.mphase), mdata(other.mdata)
    {}

    /// Update the chemical properties of the phase.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the phase (in mol)
    /// @param extra The extra properties evaluated in the activity models
    auto update(const real& T, const real& P, ArrayXrConstRef n, Map<String, Any>& extra)
    {
        _update<false>(T, P, n, extra);
    }

    /// Update the chemical properties of the phase using ideal activity models.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the phase (in mol)
    /// @param extra The extra properties evaluated in the activity models
    auto updateIdeal(const real& T, const real& P, ArrayXrConstRef n, Map<String, Any>& extra)
    {
        _update<true>(T, P, n, extra);
    }

    /// Update the chemical properties of the phase with given data.
    auto updateWithData(const ChemicalPropsPhaseBaseData<Real, Array>& data)
    {
        mdata = data;
    }

    /// Return the underlying Phase object.
    auto phase() const -> const Phase&
    {
        return mphase;
    }

    /// Return the primary chemical property data of the phase from which others are calculated.
    auto data() const -> const ChemicalPropsPhaseBaseData<Real, Array>&
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

    /// Return the amounts of the species in the phase (in mol).
    auto speciesAmounts() const -> ArrayXrConstRef
    {
        return mdata.n;
    }

    /// Return the mole fractions of the species in the phase.
    auto speciesMoleFractions() const -> ArrayXrConstRef
    {
        return mdata.x;
    }

    /// Return the ln activity coefficients of the species in the phase.
    auto speciesActivityCoefficientsLn() const -> ArrayXrConstRef
    {
        return mdata.ln_g;
    }

    /// Return the ln activities of the species in the phase.
    auto speciesActivitiesLn() const -> ArrayXrConstRef
    {
        return mdata.ln_a;
    }

    /// Return the chemical potentials of the species in the phase (in J/mol).
    auto speciesChemicalPotentials() const -> ArrayXrConstRef
    {
        return mdata.u;
    }

    /// Return the standard partial molar volumes of the species in the phase (in m³/mol).
    auto speciesStandardVolumes() const -> ArrayXrConstRef
    {
        return mdata.V0;
    }

    /// Return the temperature derivative of the standard molar volumes of the species in the phase (in m³/(mol·K)).
    auto speciesStandardVolumesT() const -> ArrayXrConstRef
    {
        return mdata.VT0;
    }

    /// Return the pressure derivative of the standard molar volumes of the species in the phase (in m³/(mol·K)).
    auto speciesStandardVolumesP() const -> ArrayXrConstRef
    {
        return mdata.VP0;
    }

    /// Return the standard partial molar Gibbs energies of formation of the species in the phase (in J/mol).
    auto speciesStandardGibbsEnergies() const -> ArrayXrConstRef
    {
        return mdata.G0;
    }

    /// Return the standard partial molar enthalpies of formation of the species in the phase (in J/mol).
    auto speciesStandardEnthalpies() const -> ArrayXrConstRef
    {
        return mdata.H0;
    }

    /// Return the standard partial molar entropies of formation of the species in the phase (in J/(mol·K)).
    auto speciesStandardEntropies() const -> ArrayXr
    {
        return (mdata.H0 - mdata.G0)/mdata.T; // from G0 = H0 - T*S0
    }

    /// Return the standard partial molar internal energies of formation of the species in the phase (in J/mol).
    auto speciesStandardInternalEnergies() const -> ArrayXr
    {
        return mdata.H0 - mdata.P * mdata.V0; // from H0 = U0 + P*V0
    }

    /// Return the standard partial molar Helmholtz energies of formation of the species in the phase (in J/mol).
    auto speciesStandardHelmholtzEnergies() const -> ArrayXr
    {
        return mdata.G0 - mdata.P * mdata.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
    }

    /// Return the standard partial molar isobaric heat capacities of the species in the phase (in J/(mol·K)).
    auto speciesStandardHeatCapacitiesConstP() const -> ArrayXrConstRef
    {
        return mdata.Cp0;
    }

    /// Return the standard partial molar isochoric heat capacities of the species in the phase (in J/(mol·K)).
    auto speciesStandardHeatCapacitiesConstV() const -> ArrayXr
    {
        return mdata.Cp0 + mdata.T * mdata.VT0 * mdata.VT0 / mdata.VP0; // from Cv0 = Cp0 + T*VT0*VT0/VP0
    }

    /// Return the molar mass of the phase (in kg/mol).
    auto molarMass() const -> real
    {
        return mdata.msum / mdata.nsum;
    }

    /// Return the molar volume of the phase (in m³/mol).
    auto molarVolume() const -> real
    {
        return (mdata.x * mdata.V0).sum() + mdata.Vx;
    }

    /// Return the temperature derivative of the molar volume of the phase (in m³/(mol·K)).
    auto molarVolumeT() const -> real
    {
        return (mdata.x * mdata.VT0).sum() + mdata.VxT;
    }

    /// Return the pressure derivative of the molar volume of the phase (in m³/(mol·Pa)).
    auto molarVolumeP() const -> real
    {
        return (mdata.x * mdata.VP0).sum() + mdata.VxP;
    }

    /// Return the molar Gibbs energy of formation of the phase (in J/mol).
    auto molarGibbsEnergy() const -> real
    {
        return (mdata.x * mdata.G0).sum() + mdata.Gx;
    }

    /// Return the molar enthalpy of formation of the phase (in J/mol).
    auto molarEnthalpy() const -> real
    {
        return (mdata.x * mdata.H0).sum() + mdata.Hx;
    }

    /// Return the molar entropy of formation of the phase (in J/(mol·K)).
    auto molarEntropy() const -> real
    {
        const auto S0 = (mdata.H0 - mdata.G0)/mdata.T; // from G0 = H0 - T*S0
        const auto Sx = (mdata.Hx - mdata.Gx)/mdata.T; // from Gx = Hx - T*Sx
        return (mdata.x * S0).sum() + Sx;
    }

    /// Return the molar internal energy of formation of the phase (in J/mol).
    auto molarInternalEnergy() const -> real
    {
        const auto U0 = mdata.H0 - mdata.P * mdata.V0; // from H0 = U0 + P*V0
        const auto Ux = mdata.Hx - mdata.P * mdata.Vx; // from Hx = U0 + P*Vx
        return (mdata.x * U0).sum() + Ux;
    }

    /// Return the molar Helmholtz energy of formation of the phase (in J/mol).
    auto molarHelmholtzEnergy() const -> real
    {
        const auto A0 = mdata.G0 - mdata.P * mdata.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
        const auto Ax = mdata.Gx - mdata.P * mdata.Vx; // from Ax = U0 - T*Sx = (Hx - P*Vx) + (Gx - Hx) = Gx - P*Vx
        return (mdata.x * A0).sum() + Ax;
    }

    /// Return the molar isobaric heat capacity of the phase (in J/(mol·K)).
    auto molarHeatCapacityConstP() const -> real
    {
        return (mdata.x * mdata.Cp0).sum() + mdata.Cpx;
    }

    /// Return the molar isochoric heat capacity of the phase (in J/(mol·K)).
    auto molarHeatCapacityConstV() const -> real
    {
        const auto Cp = molarHeatCapacityConstP();
        const auto VT = molarVolumeT();
        const auto VP = molarVolumeP();
        const auto T  = temperature();
        return Cp + T*VT*VT/VP;
    }

    /// Return the specific volume of the phase (in m³/kg).
    auto specificVolume() const -> real
    {
        return molarVolume() / molarMass();
    }

    /// Return the temperature derivative of the specific volume of the phase (in m³/(kg·K)).
    auto specificVolumeT() const -> real
    {
        return molarVolumeT() / molarMass();
    }

    /// Return the pressure derivative of the specific volume of the phase (in m³/(kg·Pa)).
    auto specificVolumeP() const -> real
    {
        return molarVolumeP() / molarMass();
    }

    /// Return the specific Gibbs energy of formation of the phase (in J/kg).
    auto specificGibbsEnergy() const -> real
    {
        return molarGibbsEnergy() / molarMass();
    }

    /// Return the specific enthalpy of formation of the phase (in J/kg).
    auto specificEnthalpy() const -> real
    {
        return molarEnthalpy() / molarMass();
    }

    /// Return the specific entropy of formation of the phase (in J/(kg·K)).
    auto specificEntropy() const -> real
    {
        return molarEntropy() / molarMass();
    }

    /// Return the specific internal energy of formation of the phase (in J/kg).
    auto specificInternalEnergy() const -> real
    {
        return molarInternalEnergy() / molarMass();
    }

    /// Return the specific Helmholtz energy of formation of the phase (in J/kg).
    auto specificHelmholtzEnergy() const -> real
    {
        return molarHelmholtzEnergy() / molarMass();
    }

    /// Return the specific isobaric heat capacity of the phase (in J/(kg·K)).
    auto specificHeatCapacityConstP() const -> real
    {
        return molarHeatCapacityConstP() / molarMass();
    }

    /// Return the specific isochoric heat capacity of the phase (in J/(kg·K)).
    auto specificHeatCapacityConstV() const -> real
    {
        return molarHeatCapacityConstV() / molarMass();
    }

    /// Return the density of the phase (in kg/m³).
    auto density() const -> real
    {
        return molarMass() / molarVolume();
    }

    /// Return the sum of species amounts in the phase (in mol).
    auto amount() const -> real
    {
        return mdata.nsum;
    }

    /// Return the sum of species masses in the phase (in kg).
    auto mass() const -> real
    {
        return mdata.msum;
    }

    /// Return the volume of the phase (in m³).
    auto volume() const -> real
    {
        return molarVolume() * amount();
    }

    /// Return the temperature derivative of the volume of the phase (in m³/K).
    auto volumeT() const -> real
    {
        return molarVolumeT() * amount();
    }

    /// Return the pressure derivative of the volume of the phase (in m³/Pa).
    auto volumeP() const -> real
    {
        return molarVolumeP() * amount();
    }

    /// Return the Gibbs energy of the phase (in J).
    auto gibbsEnergy() const -> real
    {
        return molarGibbsEnergy() * amount();
    }

    /// Return the enthalpy of the phase (in J).
    auto enthalpy() const -> real
    {
        return molarEnthalpy() * amount();
    }

    /// Return the entropy of the phase (in J/K).
    auto entropy() const -> real
    {
        return molarEntropy() * amount();
    }

    /// Return the internal energy of the phase (in J).
    auto internalEnergy() const -> real
    {
        return molarInternalEnergy() * amount();
    }

    /// Return the Helmholtz energy of the phase (in J).
    auto helmholtzEnergy() const -> real
    {
        return molarHelmholtzEnergy() * amount();
    }

    /// Return the isobaric heat capacity of the phase (in J/K).
    auto heatCapacityConstP() const -> real
    {
        return molarHeatCapacityConstP() * amount();
    }

    /// Return the isochoric heat capacity of the phase (in J/K).
    auto heatCapacityConstV() const -> real
    {
        return molarHeatCapacityConstV() * amount();
    }

    /// Return the speed of sound in the phase (in m/s).
    auto soundSpeed() const -> real
    {
        const auto rho = density();
        const auto Cp = molarHeatCapacityConstP();
        const auto Cv = molarHeatCapacityConstV();
        const auto VP = molarVolumeP();
        const auto MM = molarMass();
        return sqrt(-MM * Cp/Cv / VP) / rho;
    }

    /// Assign the given array data to this ChemicalPropsPhaseBase object.
    auto operator=(const ArrayStream<Real>& array)
    {
        mdata = array;
        return *this;
    }

    /// Convert this ChemicalPropsPhaseBase object into an array.
    operator ArrayStream<Real>() const
    {
        return mdata;
    }

    // Ensure other ChemicalPropsPhaseBase types are friend among themselves.
    template<typename RX, typename AX>
    friend class ChemicalPropsPhaseBase;

private:
    /// The phase associated with these primary chemical properties.
    Phase mphase;

    /// The primary chemical property data of the phase from which others are calculated.
    ChemicalPropsPhaseBaseData<Real, Array> mdata;

private:

    /// Update the chemical properties of the phase.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the phase (in mol)
    /// @param extra The extra data mapped to activity mode
    template<bool use_ideal_activity_model>
    auto _update(const real& T, const real& P, ArrayXrConstRef n, Map<String, Any>& extra)
    {
        mdata.T = T;
        mdata.P = P;
        mdata.n = n;

        const auto R = universalGasConstant;

        auto& nsum = mdata.nsum;
        auto& msum = mdata.msum;
        auto& x    = mdata.x;
        auto& G0   = mdata.G0;
        auto& H0   = mdata.H0;
        auto& V0   = mdata.V0;
        auto& VT0  = mdata.VT0;
        auto& VP0  = mdata.VP0;
        auto& Cp0  = mdata.Cp0;
        auto& Vx  = mdata.Vx;
        auto& VxT = mdata.VxT;
        auto& VxP = mdata.VxP;
        auto& Gx  = mdata.Gx;
        auto& Hx  = mdata.Hx;
        auto& Cpx = mdata.Cpx;
        auto& ln_g = mdata.ln_g;
        auto& ln_a = mdata.ln_a;
        auto& u    = mdata.u;

        const auto& species = phase().species();
        const auto N = species.size();

        assert(    n.size() == N );
        assert(   G0.size() == N );
        assert(   H0.size() == N );
        assert(   V0.size() == N );
        assert(  VT0.size() == N );
        assert(  VP0.size() == N );
        assert(  Cp0.size() == N );
        assert( ln_g.size() == N );
        assert( ln_a.size() == N );
        assert(    u.size() == N );

        // Compute the standard thermodynamic properties of the species in the phase.
        StandardThermoProps aux;
        for(auto i = 0; i < N; ++i)
        {
            aux = species[i].standardThermoProps(T, P);
            G0[i]  = aux.G0;
            H0[i]  = aux.H0;
            V0[i]  = aux.V0;
            VT0[i] = aux.VT0;
            VP0[i] = aux.VP0;
            Cp0[i] = aux.Cp0;
        }

        // Compute the amount of the phase
        nsum = n.sum();

        // Compute the mass of the phase
        msum = (n * mphase.speciesMolarMasses()).sum();

        // Compute the mole fractions of the species
        if(nsum == 0.0)
            x = (N == 1) ? 1.0 : 0.0;
        else x = n / nsum;

        // Ensure there are no zero mole fractions
        error(x.minCoeff() == 0.0, "Could not compute the chemical properties of phase ",
            phase().name(), " because it has one or more species with zero amounts.");

        // Compute the activity properties of the phase
        ActivityPropsRef aprops{ Vx, VxT, VxP, Gx, Hx, Cpx, ln_g, ln_a, extra };
        ActivityArgs args{ T, P, x };
        const ActivityModel& activity_model = use_ideal_activity_model ?  // IMPORTANT: Use `const ActivityModel&` here instead of `ActivityModel`, otherwise a new model is constructed without cache, and so memoization will not take effect.
            phase().idealActivityModel() : phase().activityModel();

        if(nsum == 0.0) aprops = 0.0;
        else activity_model(aprops, args);

        // Compute the chemical potentials of the species
        u = G0 + R*T*ln_a;
    }
};

/// The chemical properties of a phase and its species.
using ChemicalPropsPhase = ChemicalPropsPhaseBase<real, ArrayXr>;

/// The non-const view to the chemical properties of a phase and its species.
using ChemicalPropsPhaseRef = ChemicalPropsPhaseBase<real&, ArrayXrRef>;

/// The const view to the chemical properties of a phase and its species.
using ChemicalPropsPhaseConstRef = ChemicalPropsPhaseBase<const real&, ArrayXrConstRef>;

} // namespace Reaktoro
