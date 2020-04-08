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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declaration for ChemicalPropsPhaseDataBase
template<typename ScalarType, typename ArrayType>
struct ChemicalPropsPhaseDataBase;

/// The primary chemical property data of a phase.
using ChemicalPropsPhaseData = ChemicalPropsPhaseDataBase<real, ArrayXr>;

/// The non-const view to primary chemical property data of a phase.
using ChemicalPropsPhaseDataRef = ChemicalPropsPhaseDataBase<real&, ArrayXrRef>;

/// The const view to primary chemical property data of a phase.
using ChemicalPropsPhaseDataConstRef = ChemicalPropsPhaseDataBase<const real&, ArrayXrConstRef>;

// Forward declaration for ChemicalPropsPhaseBase
template<typename DataType>
class ChemicalPropsPhaseBase;

/// The chemical properties of a phase and its species.
using ChemicalPropsPhase = ChemicalPropsPhaseBase<ChemicalPropsPhaseData>;

/// The non-const view to the chemical properties of a phase and its species.
using ChemicalPropsPhaseRef = ChemicalPropsPhaseBase<ChemicalPropsPhaseDataRef>;

/// The const view to the chemical properties of a phase and its species.
using ChemicalPropsPhaseConstRef = ChemicalPropsPhaseBase<ChemicalPropsPhaseDataConstRef>;

/// The base type for chemical properties of a phase and its species.
template<typename DataType>
class ChemicalPropsPhaseBase
{
public:
    /// Construct a default ChemicalPropsPhaseBase instance.
    ChemicalPropsPhaseBase();

    /// Construct a ChemicalPropsPhaseBase instance.
    explicit ChemicalPropsPhaseBase(DataType data);

    /// Construct a ChemicalPropsPhaseBase instance.
    template<typename OtherDataType>
    ChemicalPropsPhaseBase(ChemicalPropsPhaseBase<OtherDataType>& props);

    /// Construct a ChemicalPropsPhaseBase instance.
    template<typename OtherDataType>
    ChemicalPropsPhaseBase(const ChemicalPropsPhaseBase<OtherDataType>& props);

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

    /// Return the volume of the phase (in m3).
    auto volume() const;

    // Ensure other ChemicalPropsPhaseBase types are friend among themselves.
    template<typename OtherDataType>
    friend class ChemicalPropsPhaseBase;

private:
    /// The primary chemical property data of the phase from which others are calculated.
    DataType props;
};

/// The base type for primary chemical property data of a phase from which others are computed.
template<typename ScalarType, typename ArrayType>
struct ChemicalPropsPhaseDataBase
{
    /// The temperature of the phase (in K).
    ScalarType T;

    /// The pressure of the phase (in Pa).
    ScalarType P;

    /// The sum of species amounts in the phase (in mol).
    ScalarType amount;

    /// The mole fractions of the species in the phase (in mol/mol).
    ArrayType x;

    /// The standard molar Gibbs energies of the species in the phase (in J/mol)
    ArrayType G0;

    /// The standard molar enthalpies of the species in the phase (in J/mol)
    ArrayType H0;

    /// The standard molar volumes of the species in the phase (in m3/mol)
    ArrayType V0;

    /// The standard molar isobaric heat capacities of the species in the phase (in J/(mol·K))
    ArrayType Cp0;

    /// The standard molar isochoric heat capacities of the species in the phase (in J/(mol·K))
    ArrayType Cv0;

    /// The excess molar volume of the phase (in m3/mol).
    ScalarType Vex;

    /// The temperature derivative of the excess molar volume at constant pressure (in m3/(mol*K)).
    ScalarType VexT;

    /// The pressure derivative of the excess molar volume at constant temperature (in m3/(mol*Pa)).
    ScalarType VexP;

    /// The excess molar Gibbs energy of the phase (in units of J/mol).
    ScalarType Gex;

    /// The excess molar enthalpy of the phase (in units of J/mol).
    ScalarType Hex;

    /// The excess molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    ScalarType Cpex;

    /// The excess molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    ScalarType Cvex;

    /// The activity coefficients (natural log) of the species in the phase.
    ArrayType ln_g;

    /// The activities (natural log) of the species in the phase.
    ArrayType ln_a;
};

template<typename DataType>
ChemicalPropsPhaseBase<DataType>::ChemicalPropsPhaseBase()
{}

template<typename DataType>
ChemicalPropsPhaseBase<DataType>::ChemicalPropsPhaseBase(DataType data)
: props(std::move(data))
{}

template<typename DataType>
template<typename OtherDataType>
ChemicalPropsPhaseBase<DataType>::ChemicalPropsPhaseBase(ChemicalPropsPhaseBase<OtherDataType>& other)
: props{
    other.props.T,
    other.props.P,
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

template<typename DataType>
template<typename OtherDataType>
ChemicalPropsPhaseBase<DataType>::ChemicalPropsPhaseBase(const ChemicalPropsPhaseBase<OtherDataType>& other)
: props{
    other.props.T,
    other.props.P,
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

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::data() const
{
    return props;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::temperature() const
{
    return props.T;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::pressure() const
{
    return props.P;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::speciesAmounts() const
{
    return props.x * props.amount;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::moleFractions() const
{
    return props.x;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::lnActivityCoefficients() const
{
    return props.ln_g;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::lnActivities() const
{
    return props.ln_a;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::chemicalPotentials() const
{
    const auto R = universalGasConstant;
    return props.G0 + R*props.T * props.ln_a;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::standardGibbsEnergies() const
{
    return props.G0;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::standardEnthalpies() const
{
    return props.H0;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::standardVolumes() const
{
    return props.V0;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::standardEntropies() const
{
    return (props.H0 - props.G0)/props.T; // from G0 = H0 - T*S0
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::standardInternalEnergies() const
{
    return props.H0 - props.P * props.V0; // from H0 = U0 + P*V0
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::standardHelmholtzEnergies() const
{
    return props.G0 - props.P * props.V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::standardHeatCapacitiesConstP() const
{
    return props.Cp0;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::standardHeatCapacitiesConstV() const
{
    return props.Cv0;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarGibbsEnergy() const
{
    return (props.x * props.G0).sum() + props.Gex;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarEnthalpy() const
{
    return (props.x * props.H0).sum() + props.Hex;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarVolume() const
{
    return (props.x * props.V0).sum() + props.Vex;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarEntropy() const
{
    const auto T = temperature();
    const auto G = molarGibbsEnergy();
    const auto H = molarEnthalpy();
    return (H - G)/T;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarInternalEnergy() const
{
    const auto P = pressure();
    const auto H = molarEnthalpy();
    const auto V = molarVolume();
    return H - P*V;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarHelmholtzEnergy() const
{
    const auto T = temperature();
    const auto U = molarInternalEnergy();
    const auto S = molarEntropy();
    return U - T*S;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarHeatCapacityConstP() const
{
    return (props.x * props.Cp0).sum() + props.Cpex;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarHeatCapacityConstV() const
{
    return (props.x * props.Cv0).sum() + props.Cvex;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::molarDensity() const
{
    return 1.0/molarVolume();
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::amount() const
{
    return props.amount;
}

template<typename DataType>
auto ChemicalPropsPhaseBase<DataType>::volume() const
{
    return molarVolume() * props.amount;
}

} // namespace Reaktoro
