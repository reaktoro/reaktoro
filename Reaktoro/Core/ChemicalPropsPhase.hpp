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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.fwd.hpp>
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

    /// The sum of species amounts in the phase (in mol).
    Real amount;

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
};

/// The base type for chemical properties of a phase and its species.
template<typename R, typename A>
class ChemicalPropsPhaseBase
{
public:
    /// Construct a ChemicalPropsPhaseBase instance.
    explicit ChemicalPropsPhaseBase(const Phase& phase);

    /// Construct a ChemicalPropsPhaseBase instance.
    ChemicalPropsPhaseBase(const Phase& phase, const ChemicalPropsPhaseBaseData<R, A>& data);

    /// Construct a ChemicalPropsPhaseBase instance.
    template<typename RX, typename AX>
    ChemicalPropsPhaseBase(ChemicalPropsPhaseBase<RX, AX>& props);

    /// Construct a ChemicalPropsPhaseBase instance.
    template<typename RX, typename AX>
    ChemicalPropsPhaseBase(const ChemicalPropsPhaseBase<RX, AX>& props);

    /// Return the primary chemical property data of the phase from which others are calculated.
    auto data() const;

    /// Return the temperature of the phase (in K).
    auto temperature() const;

    /// Return the pressure of the phase (in Pa).
    auto pressure() const;

    /// Return the amounts of the species in the phase (in mol).
    auto speciesAmounts() const;

    /// Return the mole fractions of the species in the phase.
    auto moleFractions() const;

    /// Return the ln activity coefficients of the species in the phase.
    auto lnActivityCoefficients() const;

    /// Return the ln activities of the species in the phase.
    auto lnActivities() const;

    /// Return the chemical potentials of the species (in J/mol).
    auto chemicalPotentials() const;

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

    /// Return the molar Gibbs energy of the phase (in J/mol).
    auto molarGibbsEnergy() const;

    /// Return the molar enthalpy of the phase (in J/mol).
    auto molarEnthalpy() const;

    /// Return the molar volume of the phase (in m3/mol).
    auto molarVolume() const;

    /// Return the molar entropy of the phase (in J/(mol*K)).
    auto molarEntropy() const;

    /// Return the molar internal energy of the phase (in J/mol).
    auto molarInternalEnergy() const;

    /// Return the molar Helmholtz energy of the phase (in J/mol).
    auto molarHelmholtzEnergy() const;

    /// Return the molar isobaric heat capacity of the phase (in J/(mol*K)).
    auto molarHeatCapacityConstP() const;

    /// Return the molar isochoric heat capacity of the phase (in J/(mol*K)).
    auto molarHeatCapacityConstV() const;

    /// Return the density of the phase (in kg/m3).
    auto molarDensity() const;

    /// Return the sum of species amounts in the phase (in mol).
    auto amount() const;

    /// Return the sum of species masses in the phase (in kg).
    auto mass() const;

    /// Return the volume of the phase (in m3).
    auto volume() const;

    // Ensure other ChemicalPropsPhaseBase types are friend among themselves.
    template<typename RX, typename AX>
    friend class ChemicalPropsPhaseBase;

private:
    /// The phase associated with these primary chemical properties.
    Phase phase;

    /// The primary chemical property data of the phase from which others are calculated.
    ChemicalPropsPhaseBaseData<R, A> props;
};

template<typename Real, typename Array>
ChemicalPropsPhaseBase<Real, Array>::ChemicalPropsPhaseBase(const Phase& phase)
: phase(phase)
{
    const auto numspecies = phase.species().size();

    props.n    = ArrayXr::Zero(numspecies);
    props.x    = ArrayXr::Zero(numspecies);
    props.G0   = ArrayXr::Zero(numspecies);
    props.H0   = ArrayXr::Zero(numspecies);
    props.V0   = ArrayXr::Zero(numspecies);
    props.Cp0  = ArrayXr::Zero(numspecies);
    props.Cv0  = ArrayXr::Zero(numspecies);
    props.ln_g = ArrayXr::Zero(numspecies);
    props.ln_a = ArrayXr::Zero(numspecies);
}

template<typename Real, typename Array>
ChemicalPropsPhaseBase<Real, Array>::ChemicalPropsPhaseBase(const Phase& phase, const ChemicalPropsPhaseBaseData<Real, Array>& data)
: phase(phase), props(data)
{}

template<typename Real, typename Array>
template<typename RX, typename AX>
ChemicalPropsPhaseBase<Real, Array>::ChemicalPropsPhaseBase(ChemicalPropsPhaseBase<RX, AX>& other)
: phase(other.phase), props{
    other.props.T,
    other.props.P,
    other.props.n,
    other.props.amount,
    other.props.x,
    other.props.G0,
    other.props.H0,
    other.props.V0,
    other.props.Cp0,
    other.props.Cv0,
    other.props.Vex,
    other.props.VexT,
    other.props.VexP,
    other.props.Gex,
    other.props.Hex,
    other.props.Cpex,
    other.props.Cvex,
    other.props.ln_g,
    other.props.ln_a }
{}

template<typename Real, typename Array>
template<typename RX, typename AX>
ChemicalPropsPhaseBase<Real, Array>::ChemicalPropsPhaseBase(const ChemicalPropsPhaseBase<RX, AX>& other)
: phase(other.phase), props{
    other.props.T,
    other.props.P,
    other.props.n,
    other.props.amount,
    other.props.x,
    other.props.G0,
    other.props.H0,
    other.props.V0,
    other.props.Cp0,
    other.props.Cv0,
    other.props.Vex,
    other.props.VexT,
    other.props.VexP,
    other.props.Gex,
    other.props.Hex,
    other.props.Cpex,
    other.props.Cvex,
    other.props.ln_g,
    other.props.ln_a }
{}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::data() const
{
    return props;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::temperature() const
{
    return props.T;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::pressure() const
{
    return props.P;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::speciesAmounts() const
{
    return props.n;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::moleFractions() const
{
    return props.x;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::lnActivityCoefficients() const
{
    return props.ln_g;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::lnActivities() const
{
    return props.ln_a;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::chemicalPotentials() const
{
    const auto R = universalGasConstant;
    return props.G0 + R*props.T * props.ln_a;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::standardGibbsEnergies() const
{
    return props.G0;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::standardEnthalpies() const
{
    return props.H0;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::standardVolumes() const
{
    return props.V0;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::standardEntropies() const
{
    return (props.H0 - props.G0)/props.T; // from G0 = H0 - T*S0
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::standardInternalEnergies() const
{
    return props.H0 - props.P * props.V0; // from H0 = U0 + P*V0
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::standardHelmholtzEnergies() const
{
    return props.G0 - props.P * props.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::standardHeatCapacitiesConstP() const
{
    return props.Cp0;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::standardHeatCapacitiesConstV() const
{
    return props.Cv0;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarGibbsEnergy() const
{
    return (props.x * props.G0).sum() + props.Gex;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarEnthalpy() const
{
    return (props.x * props.H0).sum() + props.Hex;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarVolume() const
{
    return (props.x * props.V0).sum() + props.Vex;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarEntropy() const
{
    const auto T = temperature();
    const auto G = molarGibbsEnergy();
    const auto H = molarEnthalpy();
    return (H - G)/T;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarInternalEnergy() const
{
    const auto P = pressure();
    const auto H = molarEnthalpy();
    const auto V = molarVolume();
    return H - P*V;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarHelmholtzEnergy() const
{
    const auto T = temperature();
    const auto U = molarInternalEnergy();
    const auto S = molarEntropy();
    return U - T*S;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarHeatCapacityConstP() const
{
    return (props.x * props.Cp0).sum() + props.Cpex;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarHeatCapacityConstV() const
{
    return (props.x * props.Cv0).sum() + props.Cvex;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::molarDensity() const
{
    return 1.0/molarVolume();
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::amount() const
{
    return props.amount;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::mass() const
{
    real sum = 0.0;
    for(auto i = 0; i < phase.species().size(); ++i)
        sum += props.n[i] * phase.species(i).molarMass();
    return sum;
}

template<typename Real, typename Array>
auto ChemicalPropsPhaseBase<Real, Array>::volume() const
{
    return molarVolume() * props.amount;
}

} // namespace Reaktoro
