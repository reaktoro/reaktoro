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
#include <Reaktoro/Common/ArrayStream.hpp>
#include <Reaktoro/Common/AutoDiff.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/Phase.hpp>

namespace Reaktoro {

/// The base type for primary chemical property data of a phase from which others are computed.
template<typename Real, typename Array>
struct PhaseChemicalPropsBaseData
{
    /// The temperature of the phase (in K).
    Real T;

    /// The pressure of the phase (in Pa).
    Real P;

    /// The amounts of each species in the phase (in mol).
    Array n;

    /// The sum of species amounts in the phase (in mol).
    Real nsum;

    /// The mole fractions of the species in the phase (in mol/mol).
    Array x;

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

    /// The excess molar volume of the phase (in m3/mol).
    Real Vex;

    /// The temperature derivative of the excess molar volume at constant pressure (in m3/(mol*K)).
    Real VexT;

    /// The pressure derivative of the excess molar volume at constant temperature (in m3/(mol*Pa)).
    Real VexP;

    /// The excess molar Gibbs energy of the phase (in units of J/mol).
    Real Gex;

    /// The excess molar enthalpy of the phase (in units of J/mol).
    Real Hex;

    /// The excess molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    Real Cpex;

    /// The excess molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    Real Cvex;

    /// The activity coefficients (natural log) of the species in the phase.
    Array ln_g;

    /// The activities (natural log) of the species in the phase.
    Array ln_a;

    /// Assign the given array data to this PhaseChemicalPropsBaseData object.
    auto operator=(const ArrayStream<Real>& array)
    {
        array.to(T, P, n, nsum, x, G0, H0, V0, Cp0, Cv0, Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a);
        return *this;
    }

    /// Convert this PhaseChemicalPropsBaseData object into an array.
    operator ArrayStream<Real>() const
    {
        return {T, P, n, nsum, x, G0, H0, V0, Cp0, Cv0, Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a};
    }
};


/// The primary chemical property data of a phase from which others are computed.
using PhaseChemicalPropsData = PhaseChemicalPropsBaseData<real, ArrayXr>;

/// The primary chemical property data of a phase from which others are computed.
using PhaseChemicalPropsDataRef = PhaseChemicalPropsBaseData<real&, ArrayXrRef>;

/// The primary chemical property data of a phase from which others are computed.
using PhaseChemicalPropsDataConstRef = PhaseChemicalPropsBaseData<const real&, ArrayXrConstRef>;


/// The base type for chemical properties of a phase and its species.
template<typename Real, typename Array>
class PhaseChemicalPropsBase
{
public:
    /// Construct a PhaseChemicalPropsBase instance.
    explicit PhaseChemicalPropsBase(const Phase& phase)
    : _phase(phase)
    {
        const auto numspecies = phase.species().size();

        _data.n    = ArrayXr::Zero(numspecies);
        _data.x    = ArrayXr::Zero(numspecies);
        _data.G0   = ArrayXr::Zero(numspecies);
        _data.H0   = ArrayXr::Zero(numspecies);
        _data.V0   = ArrayXr::Zero(numspecies);
        _data.Cp0  = ArrayXr::Zero(numspecies);
        _data.Cv0  = ArrayXr::Zero(numspecies);
        _data.ln_g = ArrayXr::Zero(numspecies);
        _data.ln_a = ArrayXr::Zero(numspecies);
    }

    /// Construct a PhaseChemicalPropsBase instance.
    PhaseChemicalPropsBase(const Phase& phase, const PhaseChemicalPropsBaseData<Real, Array>& data)
    : _phase(phase), _data(data)
    {}

    /// Construct a PhaseChemicalPropsBase instance.
    template<typename RX, typename AX>
    PhaseChemicalPropsBase(PhaseChemicalPropsBase<RX, AX>& other)
    : _phase(other.phase), _data{
        other._data.T,
        other._data.P,
        other._data.n,
        other._data.nsum,
        other._data.x,
        other._data.G0,
        other._data.H0,
        other._data.V0,
        other._data.Cp0,
        other._data.Cv0,
        other._data.Vex,
        other._data.VexT,
        other._data.VexP,
        other._data.Gex,
        other._data.Hex,
        other._data.Cpex,
        other._data.Cvex,
        other._data.ln_g,
        other._data.ln_a }
    {}

    /// Construct a PhaseChemicalPropsBase instance.
    template<typename RX, typename AX>
    PhaseChemicalPropsBase(const PhaseChemicalPropsBase<RX, AX>& other)
    : _phase(other.phase), _data{
        other._data.T,
        other._data.P,
        other._data.n,
        other._data.nsum,
        other._data.x,
        other._data.G0,
        other._data.H0,
        other._data.V0,
        other._data.Cp0,
        other._data.Cv0,
        other._data.Vex,
        other._data.VexT,
        other._data.VexP,
        other._data.Gex,
        other._data.Hex,
        other._data.Cpex,
        other._data.Cvex,
        other._data.ln_g,
        other._data.ln_a }
    {}

    /// Update the chemical properties of the phase.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the phase (in mol)
    auto update(const real& T, const real& P, ArrayXrConstRef n)
    {
        _data.T = T;
        _data.P = P;
        _data.n = n;

        auto& nsum = _data.nsum;
        auto& x    = _data.x;
        auto& G0   = _data.G0;
        auto& H0   = _data.H0;
        auto& V0   = _data.V0;
        auto& Cp0  = _data.Cp0;
        auto& Cv0  = _data.Cv0;
        auto& Vex  = _data.Vex;
        auto& VexT = _data.VexT;
        auto& VexP = _data.VexP;
        auto& Gex  = _data.Gex;
        auto& Hex  = _data.Hex;
        auto& Cpex = _data.Cpex;
        auto& Cvex = _data.Cvex;
        auto& ln_g = _data.ln_g;
        auto& ln_a = _data.ln_a;

        const auto& species = _phase.species();
        const auto size = species.size();

        assert(n.size() == size);
        assert(G0.size()   == size);
        assert(H0.size()   == size);
        assert(V0.size()   == size);
        assert(Cp0.size()  == size);
        assert(Cv0.size()  == size);
        assert(ln_g.size() == size);
        assert(ln_a.size() == size);

        // The function that evaluates the standard thermodynamic properties of the species in the phase.
        auto evalStandardThermoProps = [&]()
        {
            StandardThermoProps aux;
            for(auto i = 0; i < size; ++i)
            {
                auto standard_thermo_props_fn = _phase.species(i).standardThermoPropsFn();
                aux = standard_thermo_props_fn ? standard_thermo_props_fn(T, P) : StandardThermoProps{};
                G0[i]  = aux.G0;
                H0[i]  = aux.H0;
                V0[i]  = aux.V0;
                Cp0[i] = aux.Cp0;
                Cv0[i] = aux.Cv0;
            }
        };

        // The function that evaluates the activity properties of the phase.
        auto evalActivityProps = [&]()
        {
            nsum = n.sum();

            if(nsum == 0.0)
                x = (size == 1) ? 1.0 : 0.0;
            else x = n / nsum;

            Vec<Any> extra;
            ActivityPropsRef aprops{ Vex, VexT, VexP, Gex, Hex, Cpex, Cvex, ln_g, ln_a };
            ActivityArgs args{ T, P, x, extra };

            if(nsum == 0.0)
                aprops = 0.0;
            else _phase.activityPropsFn()(aprops, args);
        };

        evalStandardThermoProps();
        evalActivityProps();
    }

    /// Update the chemical properties of the phase.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the phase (in mol)
    /// @param wrtvar The variable with respect to automatic differentiation should be carried out.
    auto update(const real& T, const real& P, ArrayXrConstRef n, Wrt<real&> wrtvar)
    {
        autodiff::seed(wrtvar);
        update(T, P, n);
        autodiff::unseed(wrtvar);
    }

    /// Update the chemical properties of the phase.
    /// @param T The temperature condition (in K)
    /// @param P The pressure condition (in Pa)
    /// @param n The amounts of the species in the phase (in mol)
    /// @param wrtvar The index of the species amount variable for which automatic differentiation is carried out.
    auto update(const real& T, const real& P, ArrayXrRef n, Index wrtvar)
    {
        update(T, P, n, wrt(n[wrtvar]));
    }

    /// Return the underlying Phase object.
    auto phase() const -> const Phase&
    {
        return _phase;
    }

    /// Return the primary chemical property data of the phase from which others are calculated.
    auto data() const -> const PhaseChemicalPropsBaseData<Real, Array>&
    {
        return _data;
    }

    /// Return the temperature of the phase (in K).
    auto temperature() const -> real
    {
        return _data.T;
    }

    /// Return the pressure of the phase (in Pa).
    auto pressure() const -> real
    {
        return _data.P;
    }

    /// Return the amounts of the species in the phase (in mol).
    auto speciesAmounts() const -> ArrayXrConstRef
    {
        return _data.n;
    }

    /// Return the mole fractions of the species in the phase.
    auto moleFractions() const -> ArrayXrConstRef
    {
        return _data.x;
    }

    /// Return the ln activity coefficients of the species in the phase.
    auto lnActivityCoefficients() const -> ArrayXrConstRef
    {
        return _data.ln_g;
    }

    /// Return the ln activities of the species in the phase.
    auto lnActivities() const -> ArrayXrConstRef
    {
        return _data.ln_a;
    }

    /// Return the chemical potentials of the species (in J/mol).
    auto chemicalPotentials() const -> ArrayXr
    {
        const auto R = universalGasConstant;
        return _data.G0 + R*_data.T * _data.ln_a;
    }

    /// Return the standard partial molar Gibbs energies of the species (in J/mol).
    auto standardGibbsEnergies() const -> ArrayXrConstRef
    {
        return _data.G0;
    }

    /// Return the standard partial molar enthalpies of the species (in J/mol).
    auto standardEnthalpies() const -> ArrayXrConstRef
    {
        return _data.H0;
    }

    /// Return the standard partial molar volumes of the species (in m3/mol).
    auto standardVolumes() const -> ArrayXrConstRef
    {
        return _data.V0;
    }

    /// Return the standard partial molar entropies of the species (in J/(mol*K)).
    auto standardEntropies() const -> ArrayXr
    {
        return (_data.H0 - _data.G0)/_data.T; // from G0 = H0 - T*S0
    }

    /// Return the standard partial molar internal energies of the species (in J/mol).
    auto standardInternalEnergies() const -> ArrayXr
    {
        return _data.H0 - _data.P * _data.V0; // from H0 = U0 + P*V0
    }

    /// Return the standard partial molar Helmholtz energies of the species (in J/mol).
    auto standardHelmholtzEnergies() const -> ArrayXr
    {
        return _data.G0 - _data.P * _data.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
    }

    /// Return the standard partial molar isobaric heat capacities of the species (in J/(mol*K)).
    auto standardHeatCapacitiesConstP() const -> ArrayXrConstRef
    {
        return _data.Cp0;
    }

    /// Return the standard partial molar isochoric heat capacities of the species (in J/(mol*K)).
    auto standardHeatCapacitiesConstV() const -> ArrayXrConstRef
    {
        return _data.Cv0;
    }

    /// Return the molar Gibbs energy of the phase (in J/mol).
    auto molarGibbsEnergy() const -> real
    {
        return (_data.x * _data.G0).sum() + _data.Gex;
    }

    /// Return the molar enthalpy of the phase (in J/mol).
    auto molarEnthalpy() const -> real
    {
        return (_data.x * _data.H0).sum() + _data.Hex;
    }

    /// Return the molar volume of the phase (in m3/mol).
    auto molarVolume() const -> real
    {
        return (_data.x * _data.V0).sum() + _data.Vex;
    }

    /// Return the molar entropy of the phase (in J/(mol*K)).
    auto molarEntropy() const -> real
    {
        const auto T = temperature();
        const auto G = molarGibbsEnergy();
        const auto H = molarEnthalpy();
        return (H - G)/T;
    }

    /// Return the molar internal energy of the phase (in J/mol).
    auto molarInternalEnergy() const -> real
    {
        const auto P = pressure();
        const auto H = molarEnthalpy();
        const auto V = molarVolume();
        return H - P*V;
    }

    /// Return the molar Helmholtz energy of the phase (in J/mol).
    auto molarHelmholtzEnergy() const -> real
    {
        const auto T = temperature();
        const auto U = molarInternalEnergy();
        const auto S = molarEntropy();
        return U - T*S;
    }

    /// Return the molar isobaric heat capacity of the phase (in J/(mol*K)).
    auto molarHeatCapacityConstP() const -> real
    {
        return (_data.x * _data.Cp0).sum() + _data.Cpex;
    }

    /// Return the molar isochoric heat capacity of the phase (in J/(mol*K)).
    auto molarHeatCapacityConstV() const -> real
    {
        return (_data.x * _data.Cv0).sum() + _data.Cvex;
    }

    /// Return the density of the phase (in kg/m3).
    auto molarDensity() const -> real
    {
        const auto V = molarVolume();
        return V ? 1.0/V : Real(0.0);
    }

    /// Return the sum of species amounts in the phase (in mol).
    auto amount() const -> real
    {
        return _data.nsum;
    }

    /// Return the sum of species masses in the phase (in kg).
    auto mass() const -> real
    {
        Real sum = 0.0;
        for(auto i = 0; i < _phase.species().size(); ++i)
            sum += _data.n[i] * _phase.species(i).molarMass();
        return sum;
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

    /// Return the volume of the phase (in m3).
    auto volume() const -> real
    {
        return molarVolume() * amount();
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

    /// Assign the given array data to this PhaseChemicalPropsBase object.
    auto operator=(const ArrayStream<Real>& array)
    {
        _data = array;
        return *this;
    }

    /// Convert this PhaseChemicalPropsBase object into an array.
    operator ArrayStream<Real>() const
    {
        return _data;
    }

    // Ensure other PhaseChemicalPropsBase types are friend among themselves.
    template<typename RX, typename AX>
    friend class PhaseChemicalPropsBase;

private:
    /// The phase associated with these primary chemical properties.
    Phase _phase;

    /// The primary chemical property data of the phase from which others are calculated.
    PhaseChemicalPropsBaseData<Real, Array> _data;
};

/// The chemical properties of a phase and its species.
using PhaseChemicalProps = PhaseChemicalPropsBase<real, ArrayXr>;

/// The non-const view to the chemical properties of a phase and its species.
using PhaseChemicalPropsRef = PhaseChemicalPropsBase<real&, ArrayXrRef>;

/// The const view to the chemical properties of a phase and its species.
using PhaseChemicalPropsConstRef = PhaseChemicalPropsBase<const real&, ArrayXrConstRef>;

} // namespace Reaktoro
